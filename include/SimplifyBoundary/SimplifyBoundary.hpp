#ifndef SIMPLIFYBOUNDARY_HPP
#define SIMPLIFYBOUNDARY_HPP
#include "GeoJSON.hpp"
#include "stdafx.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersections.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/version.h>
#include <boost/functional/hash.hpp>

namespace PS = CGAL::Polyline_simplification_2;

namespace SimplifyBoundary {
    const double eps = 1. / (1 << 20), eps2 = eps * eps,
        eps_limit = 0.25 / (1ll << 62);

    template <typename K>
    struct SimplifyProps {
        using FT = typename K::FT;
        FT search_convex_area_thre,
            search_concave_area_thre,
            search_thre,
            del_area_thre,
            del_mdis_thre,
            del_bbs_width_thre,
            not_del_pbrate,
            theta_thre1,
            theta_thre2,
            min_circle_points_num,
            k1,
            k2,
            thre_max_multipe,
            thre_gap_node_num;
        bool isRemerge, isProcess;
        SimplifyProps() {
            search_convex_area_thre = 0.1;
            search_concave_area_thre = 0.002;
            search_thre = 0.2;
            del_area_thre = 0.025;
            del_mdis_thre = 0.05;
            del_bbs_width_thre = 0.1;
            not_del_pbrate = 0.1;
            theta_thre1 = 176;
            theta_thre2 = 160;
            min_circle_points_num = 4;
            k1 = 40;
            k2 = 0.2;
            thre_max_multipe = 4;
            thre_gap_node_num = 5;
            isRemerge = true;
            isProcess = false;
        }
    };

    template <typename K>
    struct ExpandProps {
        using FT = typename K::FT;
        std::string simplify_order;
        FT offset, tri_simplify_cost;
        ExpandProps() {
            simplify_order = "101";
            offset = 100;
            tri_simplify_cost = 0.7;
        }
        ExpandProps(const std::string& type) {
            if (type == "HOUSE")
            {
                simplify_order = "101";
                offset = 100;
                tri_simplify_cost = 0.7;
            }
            else if (type == "RAMP") {
                simplify_order = "2";
                offset = 1000;
                tri_simplify_cost = 0.7;
            }
        }
    };

    template <typename K>
    class Solver {
        using FT = typename K::FT;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Circle_2 = CGAL::Circle_2<K>;
        using Line_2 = CGAL::Line_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Stop = PS::Stop_below_count_ratio_threshold;
        using Cost = PS::Squared_distance_cost;

        struct SLNode;
        class SLList;

		using SimplifyProps = SimplifyProps<K>;
		using ExpandProps = ExpandProps<K>;

        int simplify_status;
        Polygon_with_holes_2 simplify_border(const Polygon_with_holes_2& polygon, const SimplifyProps& props);
        class SmoothMethod;
        std::pair<std::vector<Polygon_2>, std::vector<Polygon_2>> split_non_simple_polygon(const Polygon_2& polygon);
        struct point_hash {
            std::size_t operator()(const Point_2& p) const
            {
                std::size_t seed = 0;
                boost::hash_combine(seed, CGAL::to_double(p.x()));
                boost::hash_combine(seed, CGAL::to_double(p.y()));
                return seed;
            }
        };
        struct point_equal {
            bool operator()(const Point_2& p1, const Point_2& p2) const {
                return p1.x() == p2.x() && p1.y() == p2.y();
            }
        };
        struct SimplifyBoxProps;
        static FT calculateProjectionLength(const Polygon_2& polygon, const Vector_2& direction) {
            Line_2 line(Point_2(0, 0), Point_2(0, 0) + direction);

            FT minProjection = polygon[0].x() * direction.x() + polygon[0].y() * direction.y();
            FT maxProjection = minProjection;

            for (auto vertex = polygon.vertices_begin(); vertex != polygon.vertices_end(); ++vertex)
            {
                FT projection = vertex->x() * direction.x() + vertex->y() * direction.y();
                if (projection < minProjection)
                    minProjection = projection;
                if (projection > maxProjection)
                    maxProjection = projection;
            }

            return maxProjection - minProjection;
        };
        int isDelPoly(Polygon_2& poly, const Vector_2& base_v, FT& best_pbRate, FT& best_poly_area, const SimplifyBoxProps& props);
        void processPoly(std::shared_ptr<SLNode> cur, std::shared_ptr<SLNode> next, const Polygon_2& best_poly, const SimplifyBoxProps& props);
        void simplify_box(SLList& ls, const std::vector<FT>& ManhattanDistance, const SimplifyBoxProps& props);
        // Extra solver for extra cases
        inline Polygon_with_holes_2 simplify(const Polygon_with_holes_2& polygon, const SimplifyProps& Sprops, const ExpandProps& Eprops);
        inline Polygon_with_holes_2 create_exterior_offset_polygons_with_holes_2(const FT& offset, const Polygon_with_holes_2& polygon) {
            std::vector<boost::shared_ptr<Polygon_2>> out_ob_poly = CGAL::create_exterior_skeleton_and_offset_polygons_2<FT, Polygon_2, Polygon_2>(offset, polygon.outer_boundary());
            Polygon_2 new_ob = *out_ob_poly[1];
            std::cout << "new_ob orientation = " << new_ob.is_clockwise_oriented() << std::endl;
            new_ob.reverse_orientation();
            std::vector<Polygon_2> new_holes;
            for (auto it = polygon.holes_begin(); it != polygon.holes_end(); it++) {
                std::vector<boost::shared_ptr<Polygon_2>> ex_poly = CGAL::create_exterior_skeleton_and_offset_polygons_2<FT, Polygon_2, Polygon_2>(offset, *it);
                Polygon_2 temp_poly = *ex_poly[1];
                temp_poly.reverse_orientation();
                std::cout << "hole orientation = " << temp_poly.is_clockwise_oriented() << std::endl;
                new_holes.push_back(temp_poly);
            }
            Polygon_with_holes_2 new_polygon(new_ob, new_holes.begin(), new_holes.end());
            return new_polygon;
        };
        inline Polygon_with_holes_2 shrink_and_expand(const Polygon_with_holes_2& polygon, const ExpandProps& props);
        inline Polygon_with_holes_2 shrink_and_expand2(const Polygon_with_holes_2& polygon, const ExpandProps& props);
        inline Polygon_with_holes_2 tri_simplify(const Polygon_with_holes_2& polygon, const ExpandProps& props);

        struct pair_hash;
        void construct_line(std::unordered_map<int, int>& loc2degree, const int& loc, const int& last_id, std::unordered_set<std::pair<int, int>, pair_hash>& done_edges, std::vector<Point_2>& line);
        void split_to_lines(std::unordered_map<int, int>& loc2degree, std::vector<std::vector<Point_2>>& lines_group);
        struct CBSpline {
            CBSpline(void);
            ~CBSpline(void);

            static void ThreeOrderBSplineInterpolaterPt(std::vector<Point_2>& line, const int& insertNum) {
                if (line.size() <= 1 || insertNum == 0) return;

                const int len = line.size();
                std::vector<int> InsertNum(len - 1);
                for (int i = 0; i < len - 1; i++)
                    InsertNum[i] = insertNum;

                int InsertNumSum = 0;
                for (int i = 0; i < len - 1; i++)  InsertNumSum += InsertNum[i];

                std::vector<Point_2> temp(len + 2);
                for (int i = 0; i < len; i++)
                    temp[i + 1] = line[i];

                temp[0] = Point_2(2 * temp[1].x() - temp[2].x(), 2 * temp[1].y() - temp[2].y());
                temp[len + 1] = Point_2(2 * temp[len].x() - temp[len - 1].x(), 2 * temp[len].y() - temp[len - 1].y());

                Point_2 pt1, pt2, pt3, pt4;
                FT t;

                line.clear();
                line = std::vector<Point_2>(len + InsertNumSum);

                int totalnum = 0;
                for (int i = 0; i < len - 1; i++) { // insert
                    pt1 = temp[i];
                    pt2 = temp[i + 1];
                    pt3 = temp[i + 2];
                    pt4 = temp[i + 3];
                    FT dt = 1.0 / (InsertNum[i] + 1);

                    for (int j = 0; j < InsertNum[i] + 1; j++) {
                        t = dt * j;
                        FT temp_x = F03(t) * pt1.x() + F13(t) * pt2.x() + F23(t) * pt3.x() + F33(t) * pt4.x();
                        FT temp_y = F03(t) * pt1.y() + F13(t) * pt2.y() + F23(t) * pt3.y() + F33(t) * pt4.y();
                        line[totalnum] = Point_2(temp_x, temp_y);
                        totalnum++;
                    }

                    if (i == len - 2) {
                        FT temp_x = F03(1) * pt1.x() + F13(1) * pt2.x() + F23(1) * pt3.x() + F33(1) * pt4.x();
                        FT temp_y = F03(1) * pt1.y() + F13(1) * pt2.y() + F23(1) * pt3.y() + F33(1) * pt4.y();
                        line[totalnum] = Point_2(temp_x, temp_y);
                        totalnum++;
                    }
                }
            };
            static FT F03(FT t) { return 1.0 / 6 * (-t * t * t + 3 * t * t - 3 * t + 1); }
            static FT F13(FT t) { return 1.0 / 6 * (3 * t * t * t - 6 * t * t + 4); }
            static FT F23(FT t) { return 1.0 / 6 * (-3 * t * t * t + 3 * t * t + 3 * t + 1); }
            static FT F33(FT t) { return 1.0 / 6 * t * t * t; }
        };
        void simplify_line(std::vector<Point_2>& line, const SimplifyProps& props);
        void smooth_segments2(std::unordered_map<int, int>& loc2degree, std::vector<Segment_2>& result, const SimplifyProps& props, const int& insertNum);
    public:
        Solver(){}
        Solver(const std::string& geojson = "", const SimplifyProps& simplifyProps = SimplifyProps(), const ExpandProps& expandProps = ExpandProps());
		Solver(const Polygon_with_holes_2& space, const SimplifyProps& simplifyProps = SimplifyProps(), const ExpandProps& expandProps = ExpandProps());

        Polygon_with_holes_2 simplify_space;
        std::vector<Polygon_2> origin_polys;
        std::vector<Polygon_2>del_polys;
        std::vector<Point_2> seg_points;
        std::unordered_set<int> del_points;             // the points we have deleted
        static bool equal(const FT& a, const FT& b, FT epsilon = eps)
        {
            return a + epsilon > b && a < b + epsilon;
        }
        static bool less(const FT& a, const FT& b, FT epsilon = eps)
        {
            return a + epsilon <= b;
        }
        static FT Cosine(const Segment_2& a, const Segment_2& b)
        {
            auto x = a.to_vector(), y = b.to_vector();
            return (x * y) / CGAL::sqrt(x * x) / CGAL::sqrt(y * y);
        }
        static FT inner_product(const Vector_2& a, const Vector_2& b)
        {
            return a.x() * b.x() + a.y() * b.y();
        }
        static FT outer_product(const Vector_2& a, const Vector_2& b)
        {
            return a.x() * b.y() - a.y() * b.x();
        }
        static FT Cosine(const Vector_2& a, const Vector_2& b) 
        {
            return inner_product(a, b) / CGAL::approximate_sqrt(a.squared_length() * b.squared_length());
        }
        static FT Sine(const Vector_2& a, const Vector_2& b)
        {
            return outer_product(a, b) / CGAL::approximate_sqrt(a.squared_length() * b.squared_length());
        }
        static FT get_vec_len(Vector_2 vec) {
            return CGAL::approximate_sqrt(vec.squared_length());
        }
        static Vector_2 get_vec_unit(Vector_2 vec) {
            FT len = get_vec_len(vec);
            return Vector_2(vec.x() / len, vec.y() / len);
        }
        static FT get_theta_function1(const FT& a, const FT& b, const FT& x, const FT& k) {
            if (x <= 0)
                return b;
            else
                return std::max((b - a) * (std::exp(CGAL::to_double((-k) * x))) + a, (b - a) * (std::exp(CGAL::to_double(k * (x - 1)))) + a);
        };
        static FT get_cos_theta(const FT& v1_len, const FT& v2_len, const FT& a, const FT& b, const FT& k) {
            if (v1_len < 0 || v2_len < 0)
                throw "There is a error on v1_len and v2_len";
            FT x = CGAL::abs(v1_len - v2_len) / (v1_len + v2_len);
            FT cos_theta = get_theta_function1(a, b, x, k);
            //std::cout << "a = " << CGAL::to_double(a) << ", b = " << CGAL::to_double(b) << ", x = " << x << ", cos_theta = " << cos_theta << std::endl;
            return cos_theta;
        }
        static FT get_theta_function2(const FT& a, const FT& b, const FT& x, const FT& k) {
            return (b - a) * (1 - std::exp(CGAL::to_double((-k) * (x - 1)))) + a;
        }
        static FT get_cos_theta(const int& v1_num, const int& v2_num, const FT& a, const FT& b, const FT& k1, const FT& k2) {
            FT cos_theta1 = get_theta_function2(a, b, v1_num, k2);
            FT cos_theta2 = get_theta_function2(a, b, v2_num, k2);
            FT res = std::max(cos_theta1, cos_theta2);
            //std::cout << "a = " << a << ", b = " << b << ", v1_num = " << v1_num << ", v2_num = " << v2_num << ", cos_theta1 = " << cos_theta1 << ", cos_theta2 = " << cos_theta2 << ", res = " << res << std::endl;
            return res;
        }
        static FT get_cos_theta(const int& v1_num, const FT& v1_len, const int& v2_num, const FT& v2_len, const FT& a, const FT& b, const FT& k1, const FT& k2) {
            FT cos_theta1 = get_theta_function2(a, b, v1_num, k2);
            FT cos_theta2 = get_theta_function2(a, b, v2_num, k2);
            if (v1_len < 0 || v2_len < 0)
                throw "There is a error on v1_len and v2_len";
            FT x = CGAL::abs(v1_len - v2_len) / (v1_len + v2_len);
            FT cos_theta3 = get_theta_function1(a, b, x, k1);
            FT res = std::max(cos_theta3, std::max(cos_theta1, cos_theta2));    // one of Node have a lot of points and 
            //std::cout << "a = " << a << ", b = " << b << ", v1_num = " << v1_num << ", v2_num = " << v2_num << ", cos_theta1 = " << cos_theta1 << ", cos_theta2 = " << cos_theta2 << ", res = " << res << std::endl;
            return res;
        }
        static bool judge_one_side(const Line_2& temp_line, const Point_2& p1, const Point_2& p2) {
            bool is_one_side = true;
            FT s1, s2;
            if (temp_line.b() != 0) {
                s1 = p1.y() - temp_line.y_at_x(p1.x());
                s2 = p2.y() - temp_line.y_at_x(p2.x());
            }
            else {
                s1 = p1.x() - temp_line.x_at_y(p1.y());
                s2 = p2.x() - temp_line.x_at_y(p2.y());
            }
            if (s1 * s2 >= 0)
                is_one_side = true;
            else
                is_one_side = false;
            return is_one_side;
        }
        Polygon_with_holes_2 origin_space;

    private:
        /// <summary>
        /// CGAL setting for precision
        /// </summary>
        inline void CGAL_setting() const {
            CGAL::Gmpfr::set_default_precision(256);
            CGAL::set_pretty_mode(std::cout);
            CGAL::set_pretty_mode(std::cerr);
        }
    };
}

namespace SimplifyBoundary{
    template <typename K>
    Solver<K>::Solver(const std::string& geojson, const SimplifyProps& simplifyProps = SimplifyProps(), const ExpandProps& expandProps = ExpandProps()) {
        // set precision
        CGAL_setting();
        try {
            origin_space = GeoJSON::geojson_to_Pwh<K>(geojson);
            simplify_space = simplify(origin_space, simplifyProps, expandProps);
        }
        catch (const std::exception& e) {
            std::cerr << "!Error at SimplifyBoundary Module: " << e.what() << std::endl;
        }
        catch (...) {
            std::cerr << "!Caught unknown exception" << std::endl;
        }
    }

    template <typename K>
    Solver<K>::Solver(const Polygon_with_holes_2& space, const SimplifyProps& simplifyProps = SimplifyProps(), const ExpandProps& expandProps = ExpandProps()) {
        // set precision
        CGAL_setting();
        try {
            simplify_space = simplify(space, simplifyProps, expandProps);
        }
        catch (const std::exception& e) {
            std::cerr << "!Error at SimplifyBoundary Module: " << e.what() << std::endl;
        }
        catch (...) {
            std::cerr << "!Caught unknown exception" << std::endl;
        }
    }

    template <typename K>
    struct Solver<K>::SLNode
    {
        int id;
        Point_2 point;

        Vector_2 start_vec_unit, end_vec_unit; // unit vector of pre and next
        FT start_vec_len, end_vec_len;   // length of pre_vec and next_vec

        std::shared_ptr<SLNode> pre, next;
        FT avg_len;
        int line_num = 1;
        std::vector<Point_2> sh_points;
        Polygon_2 best_poly;
        FT best_pbRate = 0, best_poly_area = 0, best_cos_theta1 = -9999.0, best_cos_theta2 = -9999.0, best_mdis = 0;
        FT best_score = -99999.0;
        int best_line_sum = 0;
        std::shared_ptr<SLNode> best_prenod;

        SLNode() :pre(nullptr), next(nullptr) {}
        SLNode(int i) :id(i), pre(nullptr), next(nullptr) {}
        SLNode(int i, Point_2 p, Vector_2 vec) : id(i), point(p), pre(nullptr), next(nullptr) {
            Vector_2 vec_unit = get_vec_unit(vec);
            FT vec_len = get_vec_len(vec);
            start_vec_unit = vec_unit;
            start_vec_len = vec_len;
            end_vec_unit = vec_unit;
            end_vec_len = vec_len;
            avg_len = vec_len;
            line_num = 1;
        }
        SLNode(int i, Point_2 p, Vector_2 vec_unit, FT vec_len) :id(i), point(p), start_vec_unit(vec_unit), end_vec_len(vec_len), end_vec_unit(vec_unit), avg_len(vec_len), line_num(1), pre(nullptr), next(nullptr) {}
        void copy(const std::shared_ptr<SLNode>& nod) {
            this->id = nod->id;
            this->point = nod->point;
            this->start_vec_unit = nod->start_vec_unit;
            this->end_vec_unit = nod->end_vec_unit;
            this->start_vec_len = nod->start_vec_len;
            this->end_vec_len = nod->end_vec_len;
            this->pre = nod->pre;
            this->next = nod->next;
            this->avg_len = nod->avg_len;
            this->line_num = nod->line_num;
            this->sh_points = nod->sh_points;
            this->best_poly = nod->best_poly;
            this->best_pbRate = nod->best_pbRate;
            this->best_poly_area = nod->best_poly_area;
            this->best_cos_theta1 = nod->best_cos_theta1;
            this->best_cos_theta2 = nod->best_cos_theta2;
            this->best_mdis = nod->best_mdis;
            this->best_score = nod->best_score;
            this->best_line_sum = nod->best_line_sum;
            this->best_prenod = nod->best_prenod;
        }

        void print(const int& POINTS_SIZE, int choice = 1) {
            printf("\tNode(id = %d, point_xy = (%.6f, %.6f)", id, point.x(), point.y());
            if (pre != nullptr && (id + POINTS_SIZE - pre->id) % POINTS_SIZE > 1)
                printf(", sh_id from %d to %d", pre->id, id);
            printf("\n");
        }
        Vector_2 get_start_tangent_unit() {
            if (sh_points.size() == 0)
                return start_vec_unit;
            else if (pre->id == -1)
                return start_vec_unit;
            else {
                Point_2 p1 = pre->point, p2 = sh_points[0], p3;
                if (sh_points.size() == 1)
                    p3 = point;
                else
                    p3 = sh_points[1];
                Circle_2 curve(p1, p2, p3);
                Line_2 cp1(curve.center(), p1);
                Line_2 cp1_perp = cp1.perpendicular(p1);
                Vector_2 start_tangent_unit = get_vec_unit(cp1_perp.to_vector());
                Point_2 temp_point;
                if (cp1_perp.b() != 0)
                    temp_point = Point_2(p2.x(), cp1_perp.y_at_x(p2.x()));
                else
                    temp_point = Point_2(cp1_perp.x_at_y(p2.y()), p2.y());
                bool s1 = cp1.oriented_side(p2);
                bool s2 = cp1.oriented_side(temp_point);
                if (s1 == s2)
                    return start_tangent_unit;
                else
                    return -start_tangent_unit;
            }
        }
        Vector_2 get_end_tangent_unit() {
            int sh_size = sh_points.size();
            if (sh_size == 0)
                return end_vec_unit;
            else if (pre->id == -1)
                return end_vec_unit;
            else {
                Point_2 p1, p2 = sh_points[sh_size - 1], p3 = point;
                if (sh_size == 1)
                    p1 = pre->point;
                else
                    p1 = sh_points[sh_size - 2];
                Circle_2 curve(p1, p2, p3);
                Line_2 cp3(curve.center(), p3);
                Line_2 cp3_perp = cp3.perpendicular(p3);
                Vector_2 end_tangent_unit = get_vec_unit(cp3_perp.to_vector());
                Point_2 temp_point;
                if (cp3_perp.b() != 0)
                    temp_point = Point_2(p2.x(), cp3_perp.y_at_x(p2.x()));
                else
                    temp_point = Point_2(cp3_perp.x_at_y(p2.y()), p2.y());
                bool s1 = cp3.oriented_side(p2);
                bool s2 = cp3.oriented_side(temp_point);
                if (s1 == s2)
                    return end_tangent_unit;
                else
                    return -end_tangent_unit;
            }
        }
        void rightMerge(std::shared_ptr<SLNode>& rightNode) {
            sh_points.push_back(point);
            id = rightNode->id;
            point = rightNode->point;
            end_vec_unit = rightNode->end_vec_unit;
            end_vec_len = rightNode->end_vec_len;
            next = rightNode->next;
            avg_len = (avg_len * line_num + rightNode->avg_len * rightNode->line_num) / (line_num + rightNode->line_num);
            line_num = (line_num + rightNode->line_num);
        }
        void clearState() {
            best_pbRate = 0;
            best_poly_area = 0;
            best_cos_theta1 = -9999.0;
            best_cos_theta2 = -9999.0;
            best_mdis = 0;
            best_score = -99999.0;
            best_line_sum = 0;
            best_prenod = nullptr;
        }
    };

    template <typename K>
    class Solver<K>::SLList
    {
        std::shared_ptr<SLNode> prehead, tailnext;
        int size_count = 0;
        std::vector<Polygon_2> del_poly;

    public:
        SLList() :prehead(std::make_shared<SLNode>(-1)), tailnext(std::make_shared<SLNode>(-1)), size_count(0) {
            prehead->next = tailnext;
            tailnext->pre = prehead;
        }
        SLList(int i, Point_2 p, Vector_2 vec) : prehead(std::make_shared<SLNode>(-1)), tailnext(std::make_shared<SLNode>(-1)), size_count(1) {
            std::shared_ptr<SLNode> newNode = std::make_shared<SLNode>(i, p, vec);
            prehead->next = newNode;
            newNode->pre = prehead;
            newNode->next = tailnext;
            tailnext->pre = newNode;
        }
        SLList(const std::shared_ptr<SLNode>& newNode) : prehead(std::make_shared<SLNode>(-1)), tailnext(std::make_shared<SLNode>(-1)), size_count(1) {
            prehead->next = newNode;
            newNode->pre = prehead;
            newNode->next = tailnext;
            tailnext->pre = newNode;
        };
        SLList(const std::shared_ptr<SLNode>& headNode, const std::shared_ptr<SLNode>& tailNode) :prehead(std::make_shared<SLNode>(-1)), tailnext(std::make_shared<SLNode>(-1)) {
            std::shared_ptr<SLNode> cur = headNode;
            size_count = 1;
            while (cur != nullptr) {
                size_count++;
                cur = cur->next;
            }
            prehead->next = headNode;
            headNode->pre = prehead;
            tailNode->next = tailnext;
            tailnext->pre = tailNode;
        }
        ~SLList() {
            std::shared_ptr<SLNode> current = prehead;
            std::shared_ptr<SLNode> next;
            while (current->id != -1)
            {
                next = current->next;
                current.reset(); // Release the shared_ptr and deallocate memory
                current = next;
            }
            prehead = nullptr;
            tailnext = nullptr;
        }

        void print(const int& POINTS_SIZE, int choice = 1) {
            printf("SLList (size = %d):\n", size_count);
            std::shared_ptr<SLNode> cur = prehead->next;
            while (cur->id != -1) {
                cur->print(POINTS_SIZE);
                cur = cur->next;
            }
        }
        int size() { return size_count; }

        // # operation
        void push(std::shared_ptr<SLNode> newNode) {
            auto preNode = tailnext->pre;
            preNode->next = newNode;
            newNode->pre = preNode;
            newNode->next = tailnext;
            tailnext->pre = newNode;
            size_count++;
        }
        void pop_back() {
            if (size_count == 0)
                return;
            else if (size_count == 1) {
                prehead->next = tailnext;
                tailnext->pre = prehead;
                size_count = 0;
            }
            else {
                std::shared_ptr<SLNode> secondLastNode = tailnext->pre->pre;
                secondLastNode->next = tailnext;
                tailnext->pre = secondLastNode;
                size_count--;
            }
        }
        void pop_head() {
            if (size_count == 0)
                return;
            else if (size_count == 1) {
                prehead->next = tailnext;
                tailnext->pre = prehead;
                size_count = 0;
            }
            else {
                std::shared_ptr<SLNode> secondNode = prehead->next->next;
                prehead->next = secondNode;
                secondNode->pre = prehead;
                size_count--;
            }
        }
        std::shared_ptr<SLNode> get_tail() { return tailnext->pre; }
        std::shared_ptr<SLNode> get_head() { return prehead->next; }
        void add_del(Polygon_2& ls) {
            if (ls.size() > 2)
                del_poly.push_back(ls);
        }


        // # convert to Polygon_with_holes_2
        Polygon_2 convert() {
            Polygon_2 poly;
            poly.push_back(get_tail()->point);  // make the split polygon easier.
            std::shared_ptr<SLNode> cur = prehead->next;
            while (1) {
                for (auto p : cur->sh_points)
                    poly.push_back(p);
                if (cur->id == get_tail()->id)
                    break;
                poly.push_back(cur->point);
                cur = cur->next;
            }
            return poly;
        }
        Polygon_2 unfold() {
            Polygon_2 poly;
            std::shared_ptr<SLNode> cur = prehead->next;
            while (cur->id != -1) {
                for (auto p : cur->sh_points)
                    poly.push_back(p);
                poly.push_back(cur->point);
                cur = cur->next;
            }
            return poly;
        }
        std::vector<Polygon_2> get_del() { return del_poly; }
        void clear_del() { del_poly.clear(); }
        std::vector<Point_2> get_seg_points() {
            std::vector<Point_2> res;
            std::shared_ptr<SLNode> cur = prehead->next;
            while (cur->id != -1) {
                res.push_back(cur->point);
                cur = cur->next;
            }
            return res;
        }
        void clear_surplus_points(const FT& MAX_COS_THETA_BE_ONE_LINE) {
            std::shared_ptr<SLNode> cur = prehead->next;
            while (cur->id != -1 && cur->next->id != -1) {
                std::shared_ptr<SLNode> tempNode = cur->next;
                FT cos_theta = Cosine(-cur->end_vec_unit, tempNode->start_vec_unit);
                if (cos_theta <= MAX_COS_THETA_BE_ONE_LINE) {
                    cur->id = tempNode->id;
                    cur->sh_points.push_back(cur->point);
                    cur->point = tempNode->point;
                    cur->next = tempNode->next;
                    if (tempNode->next->id != -1)
                        tempNode->next->pre = cur;
                    cur->end_vec_unit = tempNode->end_vec_unit;
                    cur->end_vec_len = tempNode->end_vec_len;
                    continue;
                }
                cur = cur->next;
            }
        }
        void split_unnormal_circle(const int& POINTS_SIZE, const int& threshold_sh_size = 3) {
            std::shared_ptr<SLNode> cur = prehead->next;
            while (cur->id != -1) {
                if (cur->sh_points.size() != 0 && cur->sh_points.size() + 2 <= threshold_sh_size)
                    if (cur->pre->id != -1) {
                        int cnt = 0;
                        while (cnt != cur->sh_points.size()) {
                            std::shared_ptr<SLNode> preNode = cur->pre;
                            int tempid = preNode->id + 1;
                            if (tempid >= POINTS_SIZE)
                                tempid -= POINTS_SIZE;
                            std::shared_ptr<SLNode> tempNode = std::make_shared<SLNode>(tempid, cur->sh_points[cnt], Vector_2(preNode->point, cur->sh_points[cnt]));
                            preNode->next = tempNode;
                            tempNode->pre = preNode;
                            tempNode->next = cur;
                            cur->pre = tempNode;
                            cnt++;
                        }
                        cur->sh_points.clear();
                    }
                cur = cur->next;
            }
        }
        void remerge(const SimplifyBoxProps& props) {
            if (this->get_head()->id != this->get_tail()->id)
                throw "Error head node is not equal to tail node";
            std::shared_ptr<SLNode> curNode = this->get_head();  // use next is to make the first node will not be merge and then pop out.
            std::shared_ptr<SLNode> lastNode = this->get_tail()->pre;
            std::shared_ptr<SLNode> nextNode = curNode->next;
            while (curNode->next->id != -1) {
                if (size_count <= 3)    // keep each poly has at least 2 points (head is equal to tail)
                    break;
                if (curNode->sh_points.size() < 1 || nextNode->sh_points.size() < 1) {
                    lastNode = curNode;
                    curNode = nextNode;
                    nextNode = nextNode->next;
                    continue;
                }
                Vector_2 left_vec_unit = -curNode->end_vec_unit;
                Vector_2 right_vec_unit = nextNode->start_vec_unit;

                // must be one side
                Line_2 temp_line(curNode->point, left_vec_unit);
                Point_2 p1 = lastNode->point, p2 = nextNode->point;
                bool is_one_side = judge_one_side(temp_line, p1, p2);

                FT cos_theta = Cosine(left_vec_unit, right_vec_unit);
                FT cos_theta_thre = get_cos_theta(curNode->line_num, curNode->avg_len, nextNode->line_num, nextNode->avg_len, props.MAX_COS_THETA_BE_ONE_LINE, props.MAX_COS_THETA_BE_ONE_TREND, props.k1, props.k2);
                if (is_one_side == true && cos_theta <= cos_theta_thre) {
                    bool istail = nextNode->next->id == -1 ? true : false;
                    bool ishead = curNode->pre->id == -1 ? true : false;
                    curNode->sh_points.push_back(curNode->point);
                    for (auto p : nextNode->sh_points)
                        curNode->sh_points.push_back(p);
                    curNode->id = nextNode->id;
                    curNode->point = nextNode->point;
                    curNode->end_vec_unit = nextNode->end_vec_unit;
                    curNode->end_vec_len = nextNode->end_vec_len;
                    curNode->avg_len = (curNode->avg_len * curNode->line_num + nextNode->avg_len * nextNode->line_num) / (curNode->line_num + nextNode->line_num);
                    curNode->line_num += nextNode->line_num;
                    curNode->next = nextNode->next;
                    nextNode->next->pre = curNode;
                    size_count--;

                    nextNode = nextNode->next;
                    if (ishead) {
                        auto temp_pretail = this->get_tail()->pre;
                        auto temp_tailnext = this->get_tail()->next;
                        this->get_tail()->copy(this->get_head());
                        this->get_tail()->pre = temp_pretail;
                        this->get_tail()->next = temp_tailnext;
                    }
                    else if (istail) {
                        auto temp_prehead = this->get_head()->pre;
                        auto temp_headnext = this->get_head()->next;
                        this->get_head()->copy(this->get_tail());
                        this->get_head()->pre = temp_prehead;
                        this->get_head()->next = temp_headnext;
                    }

                }
                else {
                    lastNode = curNode;
                    curNode = nextNode;
                    nextNode = nextNode->next;
                }
            }
        }
    };

    template <typename K>
    inline CGAL::Polygon_with_holes_2<K> Solver<K>::simplify_border(const Polygon_with_holes_2& polygon, const SimplifyProps& props) {
        origin_polys.push_back(polygon.outer_boundary());
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
            origin_polys.push_back(*it);

        std::cout << "*****************************\tstart simplify border\t****************************" << std::endl;
        FT MIN_THETA_BE_ONE_LINE = props.theta_thre1;
        FT MIN_THETA_BE_ONE_TREND = props.theta_thre2;

        Polygon_2 outer = origin_polys[0];
        FT outer_area = outer.area();
        FT true_area = outer_area;
        for (auto hole_it = polygon.holes_begin(); hole_it != polygon.holes_end(); ++hole_it) {
            const Polygon_2& hole = *hole_it;
            true_area += hole.area();
        }
        FT multi_num = (true_area / outer_area);
        const FT max_area_threshold_of_convex_box = props.search_convex_area_thre * multi_num;
        const FT max_area_threshold_of_concave_box = props.search_concave_area_thre * multi_num;
        const FT max_area_threshold_of_del_box = props.del_area_thre * multi_num;

        FT outer_width = outer.right_vertex()->x() - outer.left_vertex()->x();
        FT outer_height = outer.top_vertex()->y() - outer.bottom_vertex()->y();
        SimplifyBoxProps SBprops;
        SBprops.MAX_SEARCH_CONVEX_AREA = true_area * max_area_threshold_of_convex_box;
        SBprops.MAX_SEARCH_CONCAVE_AREA = true_area * max_area_threshold_of_concave_box;
        SBprops.MAX_SEARCH_X = outer_width * props.search_thre;
        SBprops.MAX_SEARCH_Y = outer_height * props.search_thre;
        SBprops.MAX_DEL_AREA = true_area * max_area_threshold_of_del_box;
        SBprops.MAX_DEL_BBS_WIDTH = CGAL::approximate_sqrt(outer_width * outer_width + outer_height * outer_height) * props.del_bbs_width_thre;
        SBprops.NOT_DEL_PBRATE = props.not_del_pbrate;
        SBprops.MAX_COS_THETA_BE_ONE_LINE = std::cos(CGAL::to_double(MIN_THETA_BE_ONE_LINE * M_PI / 180));
        SBprops.MAX_COS_THETA_BE_ONE_TREND = std::cos(CGAL::to_double(MIN_THETA_BE_ONE_TREND * M_PI / 180));
        SBprops.k1 = props.k1;
        SBprops.k2 = props.k2;
        SBprops.THRESHOLD_MAX_MULTIPE = props.thre_max_multipe;
        SBprops.THRESHOLD_GAP_NODE_NUM = props.thre_gap_node_num;
        SBprops.isProcess = props.isProcess;

        printf("There are %d holes in this Polygon!\n", origin_polys.size() - 1);

        Polygon_2 outer_boundary;
        std::vector<Polygon_2>inner_holes;
        int polyid = 0;
        for (auto poly : origin_polys) {
            SBprops.POINTS_SIZE = poly.size();
            SBprops.POLYID = polyid;
            if (SBprops.POINTS_SIZE < 2)
                continue;
            printf("There are %d points in the polygon!\n", SBprops.POINTS_SIZE);

            std::vector<FT>ManhattanDistance(SBprops.POINTS_SIZE);
            ManhattanDistance[0] = 0;
            FT max_vec_len = -1, vec_len;
            int start_id = 0;
            for (int i = 1; i < SBprops.POINTS_SIZE; i++) {
                vec_len = CGAL::approximate_sqrt(CGAL::squared_distance(poly[i], poly[i - 1]));
                ManhattanDistance[i] = ManhattanDistance[i - 1] + vec_len;
                //printf("vec_len = %.2f, MD[%d] = %.2f\n", vec_len, i, ManhattanDistance[i]);
                if (max_vec_len < vec_len) {
                    max_vec_len = vec_len;
                    start_id = i;
                }
            }
            vec_len = CGAL::approximate_sqrt(CGAL::squared_distance(poly[0], poly[SBprops.POINTS_SIZE - 1]));
            ManhattanDistance[0] = ManhattanDistance[SBprops.POINTS_SIZE - 1] + vec_len;
            SBprops.MIN_DEL_MDIS = ManhattanDistance[SBprops.POINTS_SIZE - 1] * props.del_mdis_thre;
            SBprops.START_ID = start_id;
            //printf("max_vec_len = %.2f, start_id = %d\n", max_vec_len, start_id);

            // # construct the list forward
            int cnt = 1, idx = start_id, pre_idx = start_id - 1;    // start_id must be greater than 0
            Vector_2 start_vec(poly[pre_idx], poly[start_id]);
            SLList ls_forward(start_id, poly[start_id], start_vec);
            while (cnt != SBprops.POINTS_SIZE) {
                idx++;
                if (idx == SBprops.POINTS_SIZE)
                    idx = 0;
                std::shared_ptr<SLNode> lastNode = ls_forward.get_tail();
                int last_idx = lastNode->id;
                Vector_2 left_vec_unit = -lastNode->end_vec_unit;
                Vector_2 right_vec(poly[last_idx], poly[idx]);
                Vector_2 right_vec_unit = get_vec_unit(right_vec);
                FT right_vec_len = get_vec_len(right_vec);

                // must be one side
                bool is_one_side = true;
                if (lastNode->sh_points.size() > 0 && lastNode->pre->id != -1) {
                    Line_2 temp_line(lastNode->point, lastNode->end_vec_unit);
                    Point_2 p1 = lastNode->pre->point, p2 = poly[idx];
                    is_one_side = judge_one_side(temp_line, p1, p2);
                }


                std::shared_ptr<SLNode> newNode;
                // according to avg pre len, adjust threshold
                FT pre_len;
                if (pre_idx == 0)
                    pre_len = ManhattanDistance[last_idx];
                else
                    pre_len = ManhattanDistance[last_idx] - ManhattanDistance[pre_idx];
                int pre_num = last_idx - pre_idx;
                if (pre_num < 0) {
                    pre_len += ManhattanDistance[0];
                    pre_num += SBprops.POINTS_SIZE;
                }
                FT pre_avg_len = pre_len / pre_num;
                FT left_len = pre_avg_len;
                if (abs(pre_len - right_vec_len) < abs(left_len - right_vec_len))
                    left_len = pre_len;
                // TODO: according to pre max len
                //FT pre_max_len;

                // TODO: debug sample 4
                FT cos_theta = Cosine(left_vec_unit, right_vec_unit);
                FT cos_theta_thre = get_cos_theta(left_len, right_vec_len, SBprops.MAX_COS_THETA_BE_ONE_LINE, SBprops.MAX_COS_THETA_BE_ONE_TREND, SBprops.k1);
                if (pre_avg_len < 0)
                    throw "error";
                if (cos_theta <= SBprops.MAX_COS_THETA_BE_ONE_LINE || (is_one_side == true && cos_theta <= cos_theta_thre)) {
                    lastNode->end_vec_unit = right_vec_unit;
                    lastNode->end_vec_len = right_vec_len;
                    lastNode->id = idx;
                    lastNode->sh_points.push_back(lastNode->point);
                    lastNode->point = poly[idx];
                }
                else {
                    lastNode->avg_len = pre_avg_len;
                    lastNode->line_num = pre_num;
                    pre_idx = last_idx;
                    newNode = std::make_shared<SLNode>(idx, poly[idx], right_vec);
                    ls_forward.push(newNode);
                }
                cnt++;
            }
            // # construct the list backward
            // TODO: construct the list backward then merge both list
            // # merge
            SLList ls = ls_forward;
            //ls.split_unnormal_circle(POINTS_SIZE, props.min_circle_points_num);

            //ls.print(POINTS_SIZE);
            std::shared_ptr<SLNode> head_cope_node = std::make_shared<SLNode>();
            head_cope_node->copy(ls.get_head());
            ls.push(head_cope_node);

            FT temp_MAX_SEARCH_CONVEX_AREA = SBprops.MAX_SEARCH_CONVEX_AREA;
            FT temp_MAX_SEARCH_CONCAVE_AREA = SBprops.MAX_SEARCH_CONCAVE_AREA;
            // only find the concave set;
            ls.remerge(SBprops);
            SBprops.MAX_SEARCH_CONVEX_AREA = 0;
            SBprops.MAX_SEARCH_CONCAVE_AREA = temp_MAX_SEARCH_CONCAVE_AREA;
            simplify_status = 0;
            simplify_box(ls, ManhattanDistance, SBprops);
            auto dp = ls.get_del();
            del_polys.reserve(del_polys.size() + dp.size());
            del_polys.insert(del_polys.end(), dp.begin(), dp.end());
            ls.clear_del();
            // only find the convex set;
            ls.remerge(SBprops);
            SBprops.MAX_SEARCH_CONVEX_AREA = temp_MAX_SEARCH_CONVEX_AREA;
            SBprops.MAX_SEARCH_CONCAVE_AREA = 0;
            simplify_status = 1;
            simplify_box(ls, ManhattanDistance, SBprops);
            dp = ls.get_del();
            del_polys.reserve(del_polys.size() + dp.size());
            del_polys.insert(del_polys.end(), dp.begin(), dp.end());
            ls.clear_del();
            // find the nonsample polygon
            ls.remerge(SBprops);
            SBprops.MAX_SEARCH_CONVEX_AREA = temp_MAX_SEARCH_CONVEX_AREA;
            SBprops.MAX_SEARCH_CONCAVE_AREA = temp_MAX_SEARCH_CONCAVE_AREA;
            simplify_status = 2;
            simplify_box(ls, ManhattanDistance, SBprops);
            dp = ls.get_del();
            del_polys.reserve(del_polys.size() + dp.size());
            del_polys.insert(del_polys.end(), dp.begin(), dp.end());
            ls.clear_del();
            polyid++;
            ls.pop_back();
            //ls.clear_surplus_points(MAX_COS_THETA_BE_ONE_LINE);
            del_points.clear();

            // reopen list
            //ls.print(POINTS_SIZE);

            auto sp = ls.get_seg_points();
            seg_points.insert(seg_points.end(), sp.begin(), sp.end());

            Polygon_2 simplified_poly = ls.unfold();
            if (!simplified_poly.is_simple())
                throw "poly is not simple!";
            if (outer_boundary.is_empty())
                outer_boundary = simplified_poly;
            else
                inner_holes.push_back(simplified_poly);
        }
        Polygon_with_holes_2 simplified_polygon(outer_boundary, inner_holes.begin(), inner_holes.end());
        std::cout << "*****************************\tend simplify border\t*****************************" << std::endl;
        return simplified_polygon;
    }

    template <typename K>
    class Solver<K>::SmoothMethod {
        Point_2 pre_end_point, next_start_point;
        Vector_2 pre_end_vec, next_start_vec;       // start from two points
        FT pre_avg_len, next_avg_len;
        int pre_line_num, next_line_num;
        std::vector<Point_2> pre_points, next_points;    // include all pre and next points
        int MIN_ADD_LINE_NUM;
        std::shared_ptr<SLNode> curNode, preNode, nextNode;

    public:
        SmoothMethod(const std::shared_ptr<SLNode>& curNode) {
            if (curNode->pre->id == -1)
                throw "The pre node is prehead!";
            if (curNode->next->id == -1)
                throw "The next node is tailnext!";
            this->curNode = curNode;
            preNode = curNode->pre, nextNode = curNode->next;
            pre_end_point = preNode->point;
            next_start_point = curNode->point;
            pre_end_vec = preNode->end_vec_unit;
            next_start_vec = -nextNode->start_vec_unit;
            pre_points = preNode->sh_points;
            next_points = nextNode->sh_points;
            pre_points.push_back(pre_end_point);
            next_points.push_back(next_start_point);
            pre_avg_len = preNode->avg_len;
            next_avg_len = nextNode->avg_len;
            pre_line_num = preNode->line_num;
            next_line_num = nextNode->line_num;

            MIN_ADD_LINE_NUM = std::ceil(CGAL::to_double(curNode->start_vec_len / ((pre_avg_len * pre_line_num + next_avg_len * next_line_num) / (pre_line_num + next_line_num))));
        }
        std::vector<Point_2> ChaikinSolver(const FT& MAX_COS_THETA, const Point_2& outer_point = Point_2(0, 0)) {
            if (Cosine(pre_end_vec, next_start_vec) <= MAX_COS_THETA)
                return {};
            Line_2 line_a(pre_end_point, pre_end_vec);
            Line_2 line_b(next_start_point, next_start_vec);
            auto res_points = CGAL::intersection(line_a, line_b);
            Point_2 inter_point = boost::get<Point_2>(*res_points);
            Line_2 seg_line(curNode->point, curNode->end_vec_unit);
            if ((outer_point.x() == 0 && outer_point.y() == 0) || !judge_one_side(seg_line, inter_point, outer_point))
                return {};

            std::vector<Point_2> input_array = { pre_end_point, inter_point, next_start_point };
            while (input_array.size() - 1 < MIN_ADD_LINE_NUM) {
                std::vector<Point_2> output_array;
                output_array.push_back(pre_end_point);
                for (int i = 0; i < input_array.size() - 1; i++) {
                    Point_2 p0 = input_array[i];
                    Point_2 p1 = input_array[i + 1];
                    FT p0x = p0.x(), p0y = p0.y(), p1x = p1.x(), p1y = p1.y();
                    Point_2 Q(0.75 * p0x + 0.25 * p1x, 0.75 * p0y + 0.25 * p1y);
                    Point_2 R(0.25 * p0x + 0.75 * p1x, 0.25 * p0y + 0.75 * p1y);
                    output_array.push_back(Q);
                    output_array.push_back(R);
                }
                output_array.push_back(next_start_point);
                input_array = output_array;
            }
            input_array.pop_back();
            input_array.erase(input_array.begin());
            return input_array;
        }
        std::vector<Point_2> ChaikinSolver2(const FT& MAX_COS_THETA) {
            if (Cosine(pre_end_vec, next_start_vec) <= MAX_COS_THETA)
                return {};
            Line_2 line_a(pre_end_point, pre_end_vec);
            Line_2 line_b(next_start_point, next_start_vec);
            auto res_points = CGAL::intersection(line_a, line_b);
            Point_2 inter_point = boost::get<Point_2>(*res_points);
            vector<Point_2> input_array = { pre_end_point, inter_point, next_start_point };
            while (input_array.size() - 1 < MIN_ADD_LINE_NUM) {
                vector<Point_2> output_array;
                output_array.push_back(pre_end_point);
                for (int i = 0; i < input_array.size() - 1; i++) {
                    Point_2 p0 = input_array[i];
                    Point_2 p1 = input_array[i + 1];
                    FT p0x = p0.x(), p0y = p0.y(), p1x = p1.x(), p1y = p1.y();
                    Point_2 Q(0.5 * p0x + 0.5 * p1x, 0.5 * p0y + 0.5 * p1y);
                    output_array.push_back(Q);
                }
                output_array.push_back(next_start_point);
                input_array = output_array;
            }
            input_array.pop_back();
            input_array.erase(input_array.begin());
            return input_array;
        }
    };

    template <typename K>
    std::pair<std::vector<CGAL::Polygon_2<K>>, std::vector<CGAL::Polygon_2<K>>> Solver<K>::split_non_simple_polygon(const Polygon_2& polygon) {
        std::vector<Polygon_2> convex_set, concave_set;

        Segment_2 split_line((polygon.edges_end() - 1)->source(), (polygon.edges_end() - 1)->target());
        int pointNumInLine = 0;
        for (auto it = polygon.vertices_begin(); it != polygon.vertices_end(); ++it) {
            if (!split_line.has_on(*it))
                break;
            else
                pointNumInLine++;
        }
        if (pointNumInLine == polygon.size())
            return { {}, {} };


        Polygon_2 poly;
        bool is_simple = polygon.is_simple();
        for (auto ei = polygon.edges_begin(); ei != polygon.edges_end(); ++ei) {
            Segment_2 line(ei->source(), ei->target());
            poly.push_back(ei->source());
            const auto result = CGAL::intersection(split_line, line);
            if (result) {
                if (const Point_2* inter_point = boost::get<Point_2>(&*result)) {
                    if (*inter_point == ei->source()) {
                        poly.clear();
                        poly.push_back(ei->source());
                        continue;
                    }
                    else if (*inter_point == ei->target())
                        continue;
                    poly.push_back(*inter_point);
                    if (poly.is_clockwise_oriented() == false)
                        convex_set.push_back(poly);
                    else
                        concave_set.push_back(poly);
                    poly.clear();
                    poly.push_back(*inter_point);
                }
                else if (const Segment_2* inter_segment = boost::get<Segment_2>(&*result)) {
                    if (inter_segment->source() != split_line.source() && inter_segment->source() != split_line.target())
                        poly.push_back(inter_segment->source());
                    if (poly.is_clockwise_oriented() == false)
                        convex_set.push_back(poly);
                    else
                        concave_set.push_back(poly);
                    poly.clear();
                    if (inter_segment->target() != split_line.target() && inter_segment->target() != split_line.source())
                        poly.push_back(inter_segment->target());
                }
            }
        }
        bool is_simple2 = poly.is_simple();
        //poly.push_back((polygon.edges_end() - 1)->target());
        if (poly.size() > 2) {
            if (poly.is_clockwise_oriented() == false)
                convex_set.push_back(poly);
            else
                concave_set.push_back(poly);
        }
        return { convex_set, concave_set };
    }

    template <typename K>
    struct Solver<K>::SimplifyBoxProps {
        FT MAX_SEARCH_CONVEX_AREA;
        FT MAX_SEARCH_CONCAVE_AREA;
        FT MAX_SEARCH_X;
        FT MAX_SEARCH_Y;
        FT MAX_DEL_AREA;
        FT MIN_DEL_MDIS;
        FT MAX_DEL_BBS_WIDTH;
        FT NOT_DEL_PBRATE;
        FT MAX_COS_THETA_BE_ONE_LINE;
        FT MAX_COS_THETA_BE_ONE_TREND;
        FT k1;
        FT k2;
        FT THRESHOLD_MAX_MULTIPE;
        FT THRESHOLD_GAP_NODE_NUM;
        bool isProcess;
        int POLYID;
        int POINTS_SIZE;
        int START_ID;
    };

    template <typename K>
    inline int Solver<K>::isDelPoly(Polygon_2& poly, const Vector_2& base_v, FT& pbRate, FT& poly_area, const SimplifyBoxProps& props) {
        // split convex and concave set
        std::pair<std::vector<Polygon_2>, std::vector<Polygon_2>> sets = split_non_simple_polygon(poly);
        std::vector<Polygon_2> convex_set = sets.first;
        std::vector<Polygon_2> concave_set = sets.second;
        FT convex_area = 0.0, concave_area = 0.0;
        for (auto p : convex_set) {
            FT poly_area = p.area();
            convex_area += poly_area;
        }
        if (convex_set.size() > 0 && convex_area >= props.MAX_SEARCH_CONVEX_AREA)
            return 0;
        for (auto p : concave_set) {
            FT poly_area = p.area();
            concave_area += poly_area;
        }
        concave_area = -concave_area;
        if (concave_set.size() > 0 && concave_area >= props.MAX_SEARCH_CONCAVE_AREA)
            return 0;

        // Whether it intersects with other polygons
        Segment_2 split_line((poly.edges_end() - 1)->source(), (poly.edges_end() - 1)->target());
        std::unordered_set<Point_2, point_hash, point_equal> del_inter_points;
        for (auto p : del_polys) {
            for (auto ei = p.edges_begin(); ei != p.edges_end(); ++ei) {
                Segment_2 line(ei->source(), ei->target());
                const auto result = CGAL::intersection(split_line, line);
                if (result)
                    if (const Point_2* inter_point = boost::get<Point_2>(&*result))
                        del_inter_points.insert(*inter_point);
                    else if (const Segment_2* inter_segment = boost::get<Segment_2>(&*result)) {
                        del_inter_points.insert(inter_segment->source());
                        del_inter_points.insert(inter_segment->target());
                    }
            }
        }
        bool exceptSelf = (props.MAX_SEARCH_CONVEX_AREA != 0 && props.MAX_SEARCH_CONCAVE_AREA != 0);
        for (int i = 0; i < origin_polys.size(); i++) {
            if (exceptSelf && i == props.POLYID)
                continue;
            for (auto ei = origin_polys[i].edges_begin(); ei != origin_polys[i].edges_end(); ++ei) {
                Segment_2 line(ei->source(), ei->target());
                const auto result = CGAL::intersection(split_line, line);
                if (result)
                    if (const Point_2* inter_point = boost::get<Point_2>(&*result)) {
                        if (*inter_point != split_line.source() && *inter_point != split_line.target() && del_inter_points.find(*inter_point) == del_inter_points.end())
                            return 0;
                    }
                    else if (const Segment_2* inter_segment = boost::get<Segment_2>(&*result))
                        if (del_inter_points.find(inter_segment->source()) == del_inter_points.end() && del_inter_points.find(inter_segment->target()) == del_inter_points.end())
                            return 0;
            }
        }

        // the bbs is ok?
        Polygon_2 temp_poly = poly;
        Vector_2 normalized_v = base_v / CGAL::approximate_sqrt(base_v.squared_length());
        FT bbs_width = calculateProjectionLength(temp_poly, normalized_v);
        FT bbs_height = calculateProjectionLength(temp_poly, Vector_2(-normalized_v.y(), normalized_v.x()));
        FT bbs_area = bbs_width * bbs_height;
        poly_area = convex_area + concave_area;
        pbRate = poly_area / bbs_area;
        if (pbRate < props.NOT_DEL_PBRATE)
            return 0;
        return 1;
        /* */
    }

    template <typename K>
    void Solver<K>::processPoly(std::shared_ptr<SLNode> cur, std::shared_ptr<SLNode> best_next, const Polygon_2& best_poly, const SimplifyBoxProps& props) {

    }

    template <typename K>
    inline void  Solver<K>::simplify_box(SLList& ls, const std::vector<FT>& ManhattanDistance, const SimplifyBoxProps& props) {
        //printf("SLList, start id = %d\n", ls.get_head()->id);

        std::shared_ptr<SLNode> cur = ls.get_head();
        while (cur != nullptr && cur->id != -1) {
            if (cur->best_poly_area != 0 && del_points.find(cur->best_prenod->id) == del_points.end()) {
                std::shared_ptr<SLNode> best_next = cur;
                std::shared_ptr<SLNode> temp_cur = cur->best_prenod;
                FT best_poly_area = best_next->best_poly_area;
                Polygon_2 best_poly = best_next->best_poly;

                bool case1 = best_next->best_poly_area <= props.MAX_DEL_AREA;   // when the poly area is too small
                Vector_2 base_v, normalized_v, perpendicular_v;
                FT bbs_width, bbs_height;
                base_v = Vector_2(temp_cur->point, best_next->point);
                normalized_v = base_v / CGAL::approximate_sqrt(base_v.squared_length());
                perpendicular_v = Vector_2(-normalized_v.y(), normalized_v.x());
                bbs_width = calculateProjectionLength(best_poly, normalized_v);
                bbs_height = calculateProjectionLength(best_poly, perpendicular_v);

                bool case2 = best_next->best_poly_area > props.MAX_DEL_AREA && cur->best_mdis > props.MIN_DEL_MDIS && cur->best_line_sum > 20;  // when the poly may have a house structure
                bool case3 = best_next->best_poly_area > props.MAX_DEL_AREA && bbs_height / bbs_width < 1 && bbs_width < props.MAX_DEL_BBS_WIDTH;   // when the poly may be to fat

                std::shared_ptr<SLNode> origin_curnext = temp_cur->next;
                std::shared_ptr<SLNode> origin_next = best_next;
                std::shared_ptr<SLNode> nextnextNode = best_next->next;

                if (case1 || case2 || (case3 && !props.isProcess)) {
                    std::shared_ptr<SLNode> loop_cur = temp_cur->next;
                    Line_2 cur2next_perline(best_next->point, perpendicular_v);
                    Point_2 outer_point(0, 0);    // TODO: should be the outer point of bbs
                    FT max_ydis = 0.0;
                    while (loop_cur->id != best_next->id) {
                        del_points.insert(loop_cur->id);
                        FT y_dis;
                        if (cur2next_perline.b() != 0)
                            y_dis = cur2next_perline.y_at_x(loop_cur->point.x());
                        else
                            y_dis = 0;
                        if (max_ydis < y_dis) {
                            max_ydis = y_dis;
                            outer_point = loop_cur->point;
                        }
                        loop_cur = loop_cur->next;
                    }
                    if (best_poly_area < props.MAX_SEARCH_CONCAVE_AREA)     // TODO: should be other
                        outer_point = Point_2(0, 0);

                    std::shared_ptr<SLNode> cur2nextNode = std::make_shared<SLNode>(best_next->id, best_next->point, base_v);

                    temp_cur->next = cur2nextNode;
                    nextnextNode->pre = cur2nextNode;
                    cur2nextNode->next = nextnextNode;
                    cur2nextNode->pre = temp_cur;

                    SmoothMethod sm(cur2nextNode);
                    cur2nextNode->sh_points = sm.ChaikinSolver(props.MAX_COS_THETA_BE_ONE_TREND, outer_point);

                    std::shared_ptr<SLNode> new_next2curNode = std::make_shared<SLNode>(temp_cur->id, temp_cur->point, Vector_2(best_next->point, temp_cur->point));
                    new_next2curNode->sh_points = cur2nextNode->sh_points;
                    reverse(new_next2curNode->sh_points.begin(), new_next2curNode->sh_points.end());
                    origin_curnext->pre = nullptr;
                    origin_next->next = new_next2curNode;
                    new_next2curNode->pre = origin_next;
                    SLList new_box_ls(origin_curnext, new_next2curNode);
                    Polygon_2 new_del_poly = new_box_ls.convert();
                    ls.add_del(new_del_poly);

                    cur = temp_cur->next;
                }

                if (case3 && props.isProcess) {
                    std::shared_ptr<SLNode> loop_cur = temp_cur->next;
                    while (loop_cur->id != best_next->id) {
                        del_points.insert(loop_cur->id);
                        loop_cur = loop_cur->next;
                    }

                    std::shared_ptr<SLNode> cur2nextNode = std::make_shared<SLNode>(best_next->id, best_next->point, Vector_2(temp_cur->point, best_next->point));

                    temp_cur->next = cur2nextNode;
                    nextnextNode->pre = cur2nextNode;
                    cur2nextNode->next = nextnextNode;
                    cur2nextNode->pre = temp_cur;

                    SmoothMethod sm(cur2nextNode);
                    auto temp_sh_points = sm.ChaikinSolver(props.MAX_COS_THETA_BE_ONE_TREND);

                    // TODO: complete this part code
                    Point_2 CenterPoint = temp_cur->point + base_v / 2;
                    Line_2 CenterLine(CenterPoint, normalized_v);
                    FT whRate = 0.8;
                    FT process_width = bbs_height / whRate;
                    FT shrinkRate = process_width / bbs_width;
                    if (temp_sh_points.size() >= 0) {
                        temp_sh_points.clear();
                        for (auto v = best_poly.vertices_begin(); v != best_poly.vertices_end(); ++v) {
                            Point_2 p = *v;
                            Line_2 temp_line(p, perpendicular_v);
                            auto res = CGAL::intersection(temp_line, CenterLine);
                            if (Point_2* inter_point = boost::get<Point_2>(&*res)) {
                                Vector_2 shrink_vec(*inter_point, CenterPoint);
                                shrink_vec = shrink_vec * (1 - shrinkRate);
                                p = p + shrink_vec;
                                temp_sh_points.push_back(p);
                            }
                        }
                        cur2nextNode->sh_points = temp_sh_points;
                    }
                }
            }
            cur->clearState();

            //printf("id = %d\n", cur->id);
            std::shared_ptr<SLNode> next = cur->next;
            int cnt_next = 1;

            int simple_status = 0;  // 0: didn't find anyone, 1: all sets are simple, 2: all sets are nonsample
            int line_num_sum = 0;

            // loop for the most suitable next node.
            while (next != nullptr && next->id != -1 && next->next->id != -1) {
                line_num_sum += next->sh_points.size() + 1;
                std::shared_ptr<SLNode> nextnextNode = next->next;
                Vector_2 cur2next(cur->point, next->point);
                FT cur2next_len = get_vec_len(cur2next);
                if (abs(cur2next.x()) > props.MAX_SEARCH_X || abs(cur2next.y()) > props.MAX_SEARCH_Y)
                    break;

                Vector_2 cur2next_unit = get_vec_unit(cur2next);
                const FT threshold_max_multipe = props.THRESHOLD_MAX_MULTIPE;
                Vector_2 right_vec_unit = nextnextNode->start_vec_unit;
                FT right_vec_len = nextnextNode->start_vec_len;
                int num1, num2;
                if (nextnextNode->line_num == 1 && right_vec_len < (cur2next_len / threshold_max_multipe) && nextnextNode->next->id != -1 && right_vec_len < (nextnextNode->next->avg_len / threshold_max_multipe))
                    num2 = nextnextNode->next->line_num;
                else
                    num2 = nextnextNode->line_num;
                Vector_2 left_vec_unit = -cur->end_vec_unit;
                FT left_vec_len = cur->end_vec_len;
                if (cur->line_num == 1 && left_vec_len < (cur2next_len / threshold_max_multipe) && cur->pre->id != -1 && left_vec_len < (cur->pre->avg_len / threshold_max_multipe))
                    num1 = cur->pre->line_num;
                else
                    num1 = cur->line_num;

                FT Mdis, Odis;
                int headid = cur->id, tailid = next->id;
                if (headid == 0)
                    Mdis = ManhattanDistance[tailid];
                else
                    Mdis = ManhattanDistance[tailid] - ManhattanDistance[headid];

                int gap_node_num = line_num_sum - cnt_next;

                FT cos_theta1 = Cosine(left_vec_unit, cur2next_unit);
                FT cos_theta2 = Cosine(-cur2next_unit, right_vec_unit);
                FT THRESHOLD_COS_THETA = get_cos_theta(num1, num2, props.MAX_COS_THETA_BE_ONE_LINE, props.MAX_COS_THETA_BE_ONE_TREND, props.k1, props.k2);
                // cnt_next > 1 : the num of node must be greater than 1
                if (gap_node_num <= props.THRESHOLD_GAP_NODE_NUM && cnt_next > 1 && cos_theta1 <= THRESHOLD_COS_THETA && cos_theta2 <= THRESHOLD_COS_THETA) {
                    std::shared_ptr<SLNode> sh_box_head = cur->next;
                    sh_box_head->pre = nullptr;
                    std::shared_ptr<SLNode> sh_box_tail = std::make_shared<SLNode>(cur->id, cur->point, -cur2next);
                    sh_box_tail->pre = next;
                    next->next = sh_box_tail;
                    SLList box_ls(sh_box_head, sh_box_tail);
                    Polygon_2 poly = box_ls.convert();
                    bool poly_simple_status = poly.is_simple();
                    if ((simple_status == 1 && poly_simple_status == false) || poly.size() <= 2) {  // simple polygon first
                        next->next = nextnextNode;
                        sh_box_head->pre = cur;
                        break;
                    }

                    // 0: jump out(break); 1: jump next(save poly); 
                    FT pbRate, poly_area;
                    int DelStatus = isDelPoly(poly, cur2next, pbRate, poly_area, props);    // polyid represent the polygon id

                    if (DelStatus == 0) {
                        next->next = nextnextNode;
                        sh_box_head->pre = cur;
                        break;
                    }
                    else if (DelStatus == 1) {
                        if (poly_simple_status)
                            simple_status = 1;
                        else
                            simple_status = 2;
                        FT max_area = props.MAX_SEARCH_CONCAVE_AREA + props.MAX_SEARCH_CONVEX_AREA;
                        FT score = poly_area / max_area + pbRate;
                        if (next->best_score < score) {
                            next->best_score = score;
                            next->best_pbRate = pbRate;
                            next->best_poly_area = poly_area;
                            next->best_prenod = cur;
                            next->best_poly = poly;
                            next->best_cos_theta1 = cos_theta1;
                            next->best_cos_theta2 = cos_theta2;
                            next->best_mdis = Mdis;
                            next->best_line_sum = line_num_sum;
                        }

                        next->next = nextnextNode;
                        sh_box_head->pre = cur;
                    }
                }
                next = next->next;
                cnt_next++;
            }
            cur = cur->next;
        }
    }

    template <typename K>
    inline CGAL::Polygon_with_holes_2<K> Solver<K>::simplify(const Polygon_with_holes_2& polygon, const SimplifyProps& Sprops, const ExpandProps& Eprops) {
		std::cout << "*****************************\tstart simplify\t*****************************" << std::endl;
        Polygon_with_holes_2 new_polygon = polygon;
        for (char c : Eprops.simplify_order)
            switch (c) {
            case '0':
                new_polygon = simplify_border(new_polygon, Sprops);
                break;
            case '1':
                new_polygon = shrink_and_expand(new_polygon, Eprops);
                break;
            case '2':
                new_polygon = shrink_and_expand2(new_polygon, Eprops);
                break;
            case '3':
                new_polygon = tri_simplify(new_polygon, Eprops);
                break;
            }
		std::cout << "*****************************\tend simplify\t*****************************" << std::endl;
        return new_polygon;
    }

    // shrink->expand
    template <typename K>
    inline CGAL::Polygon_with_holes_2<K> Solver<K>::shrink_and_expand(const Polygon_with_holes_2& polygon, const ExpandProps& props) {
        if (props.offset <= 0) {
            std::cout << "Warning!: offset is 0." << std::endl;
            return polygon;
        }

        std::cout << "shrink_and_expand" << std::endl;
        std::cout << "offset = " << props.offset << std::endl;
        const Polygon_with_holes_2& new_polygon = polygon;
#if CGAL_VERSION_MAJOR >= 5
        /*std::cout << "CGAL version: 5.x" << std::endl;
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> in_offset_poly = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(props.offset, polygon);
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> ex_offset_poly = CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2(props.offset, *in_offset_poly[0]);
        const Polygon_with_holes_2& new_polygon = *ex_offset_poly[0];*/
#else
        std::cout << "CGAL version: 4.x" << std::endl;
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> in_offset_poly = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2<FT, Polygon_with_holes_2, Polygon_with_holes_2>(props.offset, polygon);
        const Polygon_with_holes_2& in_Polygon = *in_offset_poly[0];
        Polygon_with_holes_2 new_polygon = create_exterior_offset_polygons_with_holes_2(props.offset, in_Polygon);
#endif

        return new_polygon;
    }

    // shrink->expand->shrink
    template <typename K>
    inline CGAL::Polygon_with_holes_2<K> Solver<K>::shrink_and_expand2(const Polygon_with_holes_2& polygon, const ExpandProps& props) {
        if (props.offset <= 0) {
            std::cout << "Warning!: offset is 0." << std::endl;
            return polygon;
        }

        std::cout << "shrink_and_expand2" << std::endl;

        std::cout << "offset = " << props.offset << std::endl;
        const Polygon_with_holes_2& new_polygon = polygon;
        // cgal 4.14-3 don't support CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2 with parameter of Poly_with_holes_2
#if CGAL_VERSION_MAJOR >= 5
        /*std::vector<boost::shared_ptr<Polygon_with_holes_2>> in_offset_poly = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2<FT, Polygon_with_holes_2, Polygon_with_holes_2>(props.offset, polygon);
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> ex_offset_poly = CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2<FT, Polygon_with_holes_2, Polygon_with_holes_2>(props.offset * 2, *in_offset_poly[0]);
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> in_offset_poly2 = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2<FT, Polygon_with_holes_2, Polygon_with_holes_2>(props.offset, *ex_offset_poly[0]);
        const Polygon_with_holes_2& new_polygon = *in_offset_poly2[0];*/
#else
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> in_offset_poly = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2<FT, Polygon_with_holes_2, Polygon_with_holes_2>(props.offset, polygon);
        const Polygon_with_holes_2& in_Polygon = *in_offset_poly[0];
        Polygon_with_holes_2 ex_offset_poly = create_exterior_offset_polygons_with_holes_2(props.offset * 2, in_Polygon);
        std::vector<boost::shared_ptr<Polygon_with_holes_2>> in_offset_poly2 = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2<FT, Polygon_with_holes_2, Polygon_with_holes_2>(props.offset, ex_offset_poly);
        const Polygon_with_holes_2& new_polygon = *in_offset_poly2[0];
#endif
        std::cout << "shrink and expand done" << std::endl;
        return new_polygon;
    }

    template <typename K>
    inline CGAL::Polygon_with_holes_2<K> Solver<K>::tri_simplify(const Polygon_with_holes_2& polygon, const ExpandProps& props) {
        if (props.tri_simplify_cost >= 1) {
            std::cout << "Warning!: simplify_cost is greater than or equal to 1.";
            return polygon;
        }

        PS::Squared_distance_cost cost;

        // CGAL 4.14-3 don't support Polygon_with_holes_2 in PS:simplify
#if CGAL_VERSION_MAJOR >= 5
        Polygon_with_holes_2 new_polygon = PS::simplify(polygon, cost, Stop(CGAL::to_double(props.tri_simplify_cost)));
#else
        Polygon_2 new_ob;
        std::vector<Polygon_2> new_holes;
        new_ob = PS::simplify(polygon.outer_boundary(), cost, Stop(props.tri_simplify_cost));
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); it++)
            new_holes.push_back(PS::simplify(*it, cost, Stop(props.tri_simplify_cost)));
        Polygon_with_holes_2 new_polygon(new_ob, new_holes.begin(), new_holes.end());
#endif
        return new_polygon;
    }

    template <typename K>
    struct Solver<K>::pair_hash
    {
        template <class T1, class T2>
        size_t operator () (std::pair<T1, T2> const& pair) const
        {
            size_t h1 = hash<T1>()(pair.first);
            size_t h2 = hash<T2>()(pair.second);
            return h1 ^ h2;
        }
    };

    template <typename K>
    inline void Solver<K>::construct_line(std::unordered_map<int, int>& loc2degree, const int& loc, const int& last_loc, std::unordered_set<std::pair<int, int>, pair_hash>& done_edges, std::vector<Point_2>& line) {
        int deg = loc2degree[loc];
        line.push_back(locations[loc].point);
        if (deg == 1) {
            return;
        }
        else if (deg == 2) {
            for (std::pair<PointData*, int> nextP : locations[loc].branches) {
                if (nextP.first->is_ans == true && nextP.first->start_loc != last_loc && nextP.first->end_loc != last_loc) {
                    if (nextP.first->start_loc != loc)
                        construct_line(loc2degree, nextP.first->start_loc, loc, done_edges, line);
                    else
                        construct_line(loc2degree, nextP.first->end_loc, loc, done_edges, line);
                    return;
                }
            }
        }
        else if (deg >= 3) {
            done_edges.insert({ loc, last_loc });
            done_edges.insert({ last_loc, loc });
            return;
        }
    }

    template <typename K>
    inline void Solver<K>::split_to_lines(std::unordered_map<int, int>& loc2degree, std::vector<std::vector<Point_2>>& lines_group) {
        unordered_set<std::pair<int, int>, pair_hash> done_edges;
        for (std::pair<int, int> t : loc2degree) {
            int loc = t.first;
            int deg = t.second;
            if (deg >= 3) {
                for (std::pair<PointData*, int> edgeX : locations[loc].branches) {
                    PointData* edge = edgeX.first;
                    int next_loc = (edge->start_loc == loc) ? edge->end_loc : edge->start_loc;
                    if (edge->is_ans == false || done_edges.find({ loc, next_loc }) != done_edges.end() || done_edges.find({ next_loc, loc }) != done_edges.end()) continue;
                    std::vector<Point_2> line;
                    line.push_back(locations[loc].point);
                    construct_line(loc2degree, next_loc, loc, done_edges, line);
                    lines_group.push_back(line);
                }
            }
        }
    }

    template <typename K>
    inline void Solver<K>::simplify_line(std::vector<Point_2>& line, const SimplifyProps& props) {
        const int POINTS_SIZE = line.size();

        std::vector<FT>ManhattanDistance(POINTS_SIZE);
        ManhattanDistance[0] = 0;
        FT  vec_len;
        for (int i = 1; i < POINTS_SIZE; i++) {
            vec_len = CGAL::sqrt(CGAL::squared_distance(line[i], line[i - 1]));
            ManhattanDistance[i] = ManhattanDistance[i - 1] + vec_len;
        }

        // construct
        int start_id = 1;

        int cnt = 1, idx = start_id, pre_idx = 0;
        Vector_2 start_vec(line[pre_idx], line[idx]);
        Point_2 start_point(line[pre_idx]);
        SLList ls_forward(start_id, line[idx], start_vec);
        while (cnt != POINTS_SIZE - 1) {
            idx++;
            std::shared_ptr<SLNode> lastNode = ls_forward.get_tail();
            int last_idx = lastNode->id;
            Vector_2 left_vec_unit = -lastNode->end_vec_unit;
            Vector_2 right_vec(line[last_idx], line[idx]);
            Vector_2 right_vec_unit = get_vec_unit(right_vec);
            FT right_vec_len = get_vec_len(right_vec);

            // must be one side
            bool is_one_side = true;
            if (lastNode->sh_points.size() > 0 && lastNode->pre->id != -1) {
                Line_2 temp_line(lastNode->point, lastNode->end_vec_unit);
                Point_2 p1 = lastNode->pre->point, p2 = line[idx];
                FT s1, s2;
                if (temp_line.b() != 0) {
                    s1 = p1.y() - temp_line.y_at_x(p1.x());
                    s2 = p2.y() - temp_line.y_at_x(p2.x());
                }
                else {
                    s1 = p1.x() - temp_line.x_at_y(p1.y());
                    s2 = p1.x() - temp_line.x_at_y(p2.y());
                }
                if (s1 * s2 >= 0)
                    is_one_side = true;
                else
                    is_one_side = false;
            }

            std::shared_ptr<SLNode> newNode;
            // according to avg pre len, adjust threshold
            FT pre_len;
            if (pre_idx == 0)
                pre_len = ManhattanDistance[last_idx];
            else
                pre_len = ManhattanDistance[last_idx] - ManhattanDistance[pre_idx];
            int pre_num = last_idx - pre_idx;
            if (pre_num < 0) {
                pre_len += ManhattanDistance[0];
                pre_num += POINTS_SIZE;
            }
            FT pre_avg_len = pre_len / pre_num;
            FT left_len = pre_avg_len;
            // TODO: according to pre max len
            //FT pre_max_len;

            // TODO: debug sample 4
            FT cos_theta = Cosine(left_vec_unit, right_vec_unit);
            FT cos_theta_thre = get_cos_theta(left_len, right_vec_len, MAX_COS_THETA_BE_ONE_LINE, MAX_COS_THETA_BE_ONE_TREND, props.k1);
            if (pre_avg_len < 0)
                throw "error";
            if ((is_one_side == true && cos_theta <= cos_theta_thre)) {
                lastNode->end_vec_unit = right_vec_unit;
                lastNode->end_vec_len = right_vec_len;
                lastNode->id = idx;
                lastNode->sh_points.push_back(lastNode->point);
                lastNode->point = line[idx];
            }
            else {
                lastNode->avg_len = pre_avg_len;
                lastNode->line_num = pre_num;
                pre_idx = last_idx;
                newNode = std::make_shared<SLNode>(idx, line[idx], right_vec);
                ls_forward.push(newNode);
            }
            cnt++;
        }

        SLList ls = ls_forward;

        // # remerge the list
        std::shared_ptr<SLNode> curNode = ls.get_head();
        std::shared_ptr<SLNode> lastNode = ls.get_tail();
        std::shared_ptr<SLNode> nextNode = curNode->next;
        while (curNode->next->id != -1) {
            Vector_2 left_vec_unit = -curNode->end_vec_unit;
            Vector_2 right_vec_unit = nextNode->start_vec_unit;

            // must be one side
            bool is_one_side = true;
            Line_2 temp_line(curNode->point, right_vec_unit);
            Point_2 p1 = lastNode->point, p2 = nextNode->point;
            FT s1, s2;
            if (temp_line.b() != 0) {
                s1 = p1.y() - temp_line.y_at_x(p1.x());
                s2 = p2.y() - temp_line.y_at_x(p2.x());
            }
            else {
                s1 = p1.x() - temp_line.x_at_y(p1.y());
                s2 = p1.x() - temp_line.x_at_y(p2.y());
            }
            if (s1 * s2 >= 0)
                is_one_side = true;
            else
                is_one_side = false;

            FT cos_theta = Cosine(left_vec_unit, right_vec_unit);
            FT cos_theta_thre = get_cos_theta(curNode->line_num, nextNode->line_num, MAX_COS_THETA_BE_ONE_LINE, MAX_COS_THETA_BE_ONE_TREND, props.k1, props.k2);
            if (is_one_side == true && cos_theta <= cos_theta_thre) {
                curNode->sh_points.push_back(curNode->point);
                for (auto p : nextNode->sh_points)
                    curNode->sh_points.push_back(p);
                curNode->id = nextNode->id;
                curNode->point = nextNode->point;
                curNode->end_vec_unit = nextNode->end_vec_unit;
                curNode->end_vec_len = nextNode->end_vec_len;
                curNode->avg_len = (curNode->avg_len * curNode->line_num + nextNode->avg_len * nextNode->line_num) / (curNode->line_num + nextNode->line_num);
                curNode->line_num += nextNode->line_num;
                curNode->next = nextNode->next;
                nextNode->next->pre = curNode;

                nextNode = nextNode->next;
            }
            else {
                lastNode = curNode;
                curNode = nextNode;
                nextNode = nextNode->next;
            }
            if (nextNode->id == -1)
                nextNode = ls.get_head();
        }
        ls.split_unnormal_circle(POINTS_SIZE, props.min_circle_points_num);

        // # simplify
        std::shared_ptr<SLNode> cur = ls.get_head();
        while (cur->id != -1) {
            std::shared_ptr<SLNode> next = cur->next;
            int cnt_next = 1;
            FT best_cos_theta1, best_cos_theta2;
            std::shared_ptr<SLNode> best_next = cur->next;
            Vector_2 temp_cur2next;

            while (next->id != -1 && next->next->id != -1) {
                std::shared_ptr<SLNode> nextnextNode = next->next;
                Vector_2 cur2next(cur->point, next->point);
                FT cur2next_len = get_vec_len(cur2next);
                if (cur2next_len > MAX_SEARCH_RADIUS)
                    break;

                Vector_2 cur2next_unit = get_vec_unit(cur2next);
                const FT threshold_max_multipe = props.thre_max_multipe;
                Vector_2 right_vec_unit = nextnextNode->start_vec_unit;
                FT right_vec_len = nextnextNode->start_vec_len;
                FT num1, num2;
                if (nextnextNode->line_num == 1 && right_vec_len < (cur2next_len / threshold_max_multipe) && nextnextNode->next->id != -1 && right_vec_len < (nextnextNode->next->avg_len / threshold_max_multipe))
                    num2 = nextnextNode->next->line_num;
                else
                    num2 = nextnextNode->line_num;
                Vector_2 left_vec_unit = -cur->end_vec_unit;
                FT left_vec_len = cur->end_vec_len;
                if (cur->line_num == 1 && left_vec_len < (cur2next_len / threshold_max_multipe) && cur->pre->id != -1 && left_vec_len < (cur->pre->avg_len / threshold_max_multipe))
                    num1 = cur->pre->line_num;
                else
                    num1 = cur->line_num;

                FT Mdis, Odis;
                int headid = cur->id, tailid = next->id;
                if (headid == 0)
                    Mdis = ManhattanDistance[tailid];
                else
                    Mdis = ManhattanDistance[tailid] - ManhattanDistance[headid];
                int line_num_gap = tailid - headid;
                if (line_num_gap < 0) {
                    Mdis += ManhattanDistance[0];
                    line_num_gap += POINTS_SIZE;
                }
                Odis = cur2next_len;
                int gap_node_num = line_num_gap - cnt_next;

                FT cos_theta1 = Cosine(left_vec_unit, cur2next_unit);
                FT cos_theta2 = Cosine(-cur2next_unit, right_vec_unit);
                FT THRESHOLD_COS_THETA = get_cos_theta(num1, num2, MAX_COS_THETA_BE_ONE_LINE, MAX_COS_THETA_BE_ONE_TREND, props.k1, props.k2);
                if (gap_node_num < props.thre_gap_node_num && cos_theta1 <= THRESHOLD_COS_THETA && cos_theta2 <= THRESHOLD_COS_THETA) {
                    best_next = next;
                    temp_cur2next = cur2next;
                    best_cos_theta1 = cos_theta1;
                    best_cos_theta2 = cos_theta2;
                }
                next = next->next;
                cnt_next++;
            }

            if (cur->next->id != best_next->id) {
                std::shared_ptr<SLNode> nextnextNode = best_next->next;
                std::shared_ptr<SLNode> tempNode = std::make_shared<SLNode>(best_next->id, best_next->point, temp_cur2next);
                cur->next = tempNode;
                nextnextNode->pre = tempNode;
                tempNode->next = nextnextNode;
                tempNode->pre = cur;

                SmoothMethod sm(tempNode);
                tempNode->sh_points = sm.ChaikinSolver(MAX_COS_THETA_BE_ONE_TREND);

                if (best_cos_theta2 <= MAX_COS_THETA_BE_ONE_LINE)
                    cur = cur->next;
            }
            cur = cur->next;
        }

        // # convert
        line.clear();
        std::shared_ptr<SLNode> tempcur = ls.get_head();
        line.push_back(start_point);
        while (tempcur->id != -1) {
            for (auto p : tempcur->sh_points)
                line.push_back(p);
            line.push_back(tempcur->point);
            tempcur = tempcur->next;
        }
    }

    template <typename K>
    inline void Solver<K>::smooth_segments2(std::unordered_map<int, int>& loc2degree, std::vector<Segment_2>& result, const SimplifyProps& props, const int& insertNum) {
        std::vector<std::vector<Point_2>> lines_group;
        split_to_lines(loc2degree, lines_group);

        for (int i = 0; i < lines_group.size(); i++) {
            std::vector<Point_2> line = lines_group[i];

            //simplify_line(line, props);
            CBSpline::ThreeOrderBSplineInterpolaterPt(line, insertNum);

            int len = line.size();
            for (int i = 1; i < len; i++) {
                result.emplace_back(line[i - 1], line[i]);
            }
        }
    }
}
#endif // SIMPLIFYBOUNDARY_HPP