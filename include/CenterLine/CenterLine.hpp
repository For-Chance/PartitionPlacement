#ifndef CENTERLINE_HPP
#define CENTERLINE_HPP
#include "GeoJSON.hpp"
#include "stdafx.h"

#include "convertKernel.hpp"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Straight_skeleton_builder_2.h>

namespace CenterLine {
    const double eps = 1. / (1 << 20), eps2 = eps * eps,
        eps_limit = 0.25 / (1ll << 62);

    struct PairHash {
        std::size_t operator()(const std::pair<int, int>& p) const
        {
            return std::hash<std::uint64_t>{}(std::uint64_t(p.first) << 32 | p.second);
        }
    };

    template <typename K>
    struct CenterLineProps {
        using FT = typename K::FT;
        FT theta_thre3;
        bool isSmooth;
        int insertNum;
        FT min_convex_height;

        CenterLineProps() {
            theta_thre3 = 150;
            isSmooth = false;
            insertNum = 10;
            min_convex_height = 0.05;
        }
    };

    using IK = CGAL::Epeck;
    using IFT = IK::FT;
    inline void get_intersection(const CGAL::Polygon_2<IK>& poly, const CGAL::Ray_2<IK>& ray, IK::Point_2& new_p, IFT& cur_dis) {
        CGAL::Point_2<IK>* tmp_p;
        CGAL::Segment_2<IK>* tmp_seg;
        for (auto it = poly.edges_begin(); it != poly.edges_end(); ++it) {
            auto inter = CGAL::intersection(*it, ray);
            if (!inter) continue;
            if ((tmp_p = boost::get<CGAL::Point_2<IK>>(&*inter))) {
                IFT dis = CGAL::squared_distance(*tmp_p, ray.source());
                if (cur_dis < 0 || dis < cur_dis) {
                    cur_dis = dis;
                    new_p = *tmp_p;
                }
            }
            else if ((tmp_seg = boost::get<CGAL::Segment_2<IK>>(&*inter))) {
                IFT dis_0 = CGAL::squared_distance(tmp_seg->source(), ray.source());
                IFT dis_1 = CGAL::squared_distance(tmp_seg->target(), ray.source());
                if (cur_dis < 0 || CGAL::min(dis_0, dis_1) < cur_dis) {
                    if (dis_0 < dis_1) {
                        cur_dis = dis_0;
                        new_p = tmp_seg->source();
                    }
                    else {
                        cur_dis = dis_1;
                        new_p = tmp_seg->target();
                    }
                }
            }
        }
    }
    inline void get_intersection2(const CGAL::Polygon_2<IK>& poly, const CGAL::Ray_2<IK>& ray, IK::Point_2& new_p, IFT& cur_dis) {
        CGAL::Point_2<IK>* tmp_p;
        CGAL::Segment_2<IK>* tmp_seg;
        for (auto it = poly.edges_begin(); it != poly.edges_end(); ++it) {
            auto inter = CGAL::intersection(*it, ray);
            if (!inter) continue;
            if ((tmp_p = boost::get<CGAL::Point_2<IK>>(&*inter))) {
                IFT dis = CGAL::squared_distance(*tmp_p, ray.source());
                if (cur_dis < 0 || dis < cur_dis) {
                    cur_dis = dis;
                    new_p = *tmp_p;
                }
            }
            else if ((tmp_seg = boost::get<CGAL::Segment_2<IK>>(&*inter))) {
                IFT dis_0 = CGAL::squared_distance(tmp_seg->source(), ray.source());
                IFT dis_1 = CGAL::squared_distance(tmp_seg->target(), ray.source());
                if (cur_dis < 0 || CGAL::min(dis_0, dis_1) < cur_dis) {
                    if (dis_0 < dis_1) {
                        cur_dis = dis_0;
                        new_p = tmp_seg->source();
                    }
                    else {
                        cur_dis = dis_1;
                        new_p = tmp_seg->target();
                    }
                }
            }
        }
    }
    template <typename K>
    inline typename K::Point_2 get_intersection(const CGAL::Polygon_with_holes_2<K>& poly, typename K::Point_2 point, typename K::Vector_2 dir) {
        KernelConverter::KernelConverter<K, CGAL::Epeck, KernelConverter::NumberConverter<K::FT, IK::FT>> to_exact;
        KernelConverter::KernelConverter<CGAL::Epeck, K, KernelConverter::NumberConverter<CGAL::Epeck::FT, K::FT, 256>> to_Gmpfr;
        CGAL::Polygon_with_holes_2<IK> space = to_exact.convert(poly);
        CGAL::Point_2<IK> p = to_exact(point), new_p;
        IFT cur_dis = -1;
        CGAL::Vector_2<IK> v = to_exact(dir);
        CGAL::Ray_2<IK> ray(p, v);
        std::cout << "polygon = " << space << std::endl << "ray = " << ray << std::endl;
        get_intersection(space.outer_boundary(), ray, new_p, cur_dis);
        for (auto h_it = space.holes_begin(); h_it != space.holes_end(); ++h_it)
            get_intersection(*h_it, ray, new_p, cur_dis);
        if (cur_dis < 0) {
            std::cerr << "get_intersection(poly, point, dir): no intersection" << std::endl << point << " " << dir << std::endl;
            throw("get_intersection(poly, point, dir): no intersection");
        }
        return to_Gmpfr(new_p);
    }
    template <typename K>
    inline typename K::Point_2 get_intersection2(const CGAL::Polygon_with_holes_2<K>& poly, typename K::Point_2 point, typename K::Vector_2 dir) {
        KernelConverter::KernelConverter<K, CGAL::Epeck, KernelConverter::NumberConverter<K::FT, IK::FT>> to_exact;
        KernelConverter::KernelConverter<CGAL::Epeck, K, KernelConverter::NumberConverter<CGAL::Epeck::FT, K::FT, 256>> to_Gmpfr;
        CGAL::Polygon_with_holes_2<IK> space = to_exact.convert(poly);
        CGAL::Point_2<IK> p = to_exact(point), new_p;
        IFT cur_dis = -1;
        CGAL::Vector_2<IK> v = to_exact(dir);
        CGAL::Ray_2<IK> ray(p, v);
        //std::cout << "polygon = " << space << std::endl << "ray = " << ray << std::endl;
        get_intersection2(space.outer_boundary(), ray, new_p, cur_dis);
        for (auto h_it = space.holes_begin(); h_it != space.holes_end(); ++h_it)
            get_intersection2(*h_it, ray, new_p, cur_dis);
        if (cur_dis < 0) {
            std::cerr << "get_intersection(poly, point, dir): no intersection" << std::endl << point << " " << dir << std::endl;
            throw("get_intersection(poly, point, dir): no intersection");
        }
        return to_Gmpfr(new_p);
    }

    template<typename K>
    inline bool evalDirection(CGAL::Vector_2<K> from, CGAL::Vector_2<K> base_v, CGAL::Vector_2<K>& used_vector, typename K::FT& value) {
        using FT = typename K::FT;
        using Solver = CenterLineSolver<K, CGAL::Polygon_with_holes_2<K>, CGAL::Polygon_2<K>>;
        bool upd = false;
        FT Cos = Solver::Cosine(base_v, from);
        if (Solver::equal(Cos, 1) || Solver::equal(Cos, 0)) {
            if (value < 1) {
                used_vector = base_v;
                value = 1;
                upd = true;
            }
        }
        else if (Cos >= 0) {
            FT t0 = 1 - Cos * Cos;
            FT t1 = Cos * Cos / 2;
            if (t0 > t1) {
                if (value < t0) {
                    used_vector = from;
                    value = t0;
                    upd = true;
                }
            }
            else {
                if (value < t1) {
                    used_vector = base_v;
                    value = t1;
                    upd = true;
                }
            }
        }
        else if (value < Cos) {
            used_vector = base_v;
            value = Cos;
            upd = true;
        }
        return upd;
    }
    template<typename K>
    inline typename K::FT evalDirection(CGAL::Vector_2<K> from, CGAL::Vector_2<K> base_v, CGAL::Vector_2<K> to) {

    }

    template <typename K>
    class Solver {
        using FT = typename K::FT;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Ray_2 = CGAL::Ray_2<K>;
        using Segment_2 = CGAL::Segment_2<K>;
        using Ss = CGAL::Straight_skeleton_2<K>;
        using SsBuilderTraits = CGAL::Straight_skeleton_builder_traits_2<K>;
        using SsBuilder = CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss>;
        using Vertex_handle = typename Ss::Vertex_handle;
        using Vertex_iterator = typename Ss::Vertex_iterator;
        using Halfedge_iterator = typename Ss::Halfedge_iterator;
        using Halfedge_handle = typename Ss::Halfedge_handle;

        using CenterLineProps = CenterLineProps<K>;
    public:
        Polygon_with_holes_2 origin_space;

        size_t seg_cnt_before_connect;
        Solver(const std::string& geojson = "", const CenterLineProps& centerlineProps = CenterLineProps());
        const std::vector<Segment_2>& centerline() const { return res_segments; }
        const std::vector<Segment_2>& sub_centerline() const { return sub_segments; }
        const std::vector<Segment_2>& smooth_centerline() const { return smooth_segments; }
        std::string centerline_geojson() const { return GeoJSON::segments_to_geojson<K>(res_segments); }
        std::string sub_centerline_geojson() const { return GeoJSON::segments_to_geojson<K>(sub_segments); }
        std::string centerline_smooth_geojson() const { return GeoJSON::segments_to_geojson<K>(smooth_segments); }
        std::vector<Point_2> get_log_points() const { return log_points; }

    private:
        /// <summary>
        /// CGAL setting for precision
        /// </summary>
        inline void CGAL_setting() const {
            CGAL::Gmpfr::set_default_precision(256);
            CGAL::set_pretty_mode(std::cout);
            CGAL::set_pretty_mode(std::cerr);
        }

        struct Location;
        struct PointData;

        struct Location {
            Point_2 point;
            FT time;
            std::vector<std::pair<PointData*, bool>> branches;
            Location() {}
            Location(Point_2 p, FT t) : point(p), time(t) {}
            int get_ans_num() {
                int ans_num = 0;
                for (auto b : branches)
                    if (b.first->is_ans == true)
                        ans_num++;
                return ans_num;
            }
            PointData* get_next_ans_edge() {
                if (get_ans_num() != 1)
                    return nullptr;
                for (auto b : branches)
                    if (b.first->is_ans == true)
                        return b.first;
                return nullptr;
            }
        };
        struct PointData {
            std::vector<Location>& locations;
            int start_loc, end_loc; // start time & location; end time & location;
            int point_id;
            Vector_2 src_vector, dest_vector; // src_vector; dest_vector; (unit vector)
            //PointData *branches[2]; // branches[2]: PointData. Points created when splitting
            std::vector<std::pair<PointData*, bool>> branches[2]; // 0: start; 1: end
            FT squared_cos_theta, inner;
            PointData* prev, * next;
            bool is_ans;

            PointData(std::vector<Location>& loc, int _start_loc,
                      const Vector_2& src, const Vector_2& dest)
                : locations(loc), start_loc(_start_loc), end_loc(-1), src_vector(src), dest_vector(dest), prev(this), next(this)
            {
                std::cout << "start_loc = " << _start_loc << std::endl;
                locations[_start_loc].branches.emplace_back(this, 0);
                FT tmp_inner_product = inner_product(src_vector, dest_vector);
                FT absolute_scale = src_vector.squared_length() * dest_vector.squared_length();
                is_ans = (outer_product(src_vector, dest_vector) >= 0 && tmp_inner_product < 0 &&
                          tmp_inner_product * tmp_inner_product > FT(0.25 - eps) * absolute_scale);
            }
            PointData(std::vector<Location>& loc, int _start_loc, int _end_loc, int id)
                : locations(loc), start_loc(_start_loc), end_loc(_end_loc), is_ans(true), prev(this), next(this), point_id(id) {}
            void set_end(int loc)
            {
                end_loc = loc;
                locations[loc].branches.emplace_back(this, 1);
            }
            Point_2 source() const { return locations[start_loc].point; }
            Point_2 target() const
            {
                if (end_loc != -1)
                    return locations[end_loc].point;
                else
                    return Point_2(FT(0) / FT(0), FT(0) / FT(0));
            }
            Point_2 point(bool port) { return port ? target() : source(); }
            int location(bool port) { return port ? end_loc : start_loc; }
            FT start_time() const { return locations[start_loc].time; }
            FT end_time() const {
                if (end_loc != -1) return locations[end_loc].time;
                else return FT(1) / FT(0);
            }
            Vector_2 speed() const
            {
                // TODO: 处理极端情况(夹角接近+-180°)
                return calc_speed(src_vector, dest_vector);
            }
            Ray_2 path() const
            {
                return Ray_2(source(), speed());
            }
            Point_2 new_loc(FT time) const
            {
                return source() + calc_real_speed(src_vector, dest_vector) * (CGAL::sqrt(time) - CGAL::sqrt(start_time()));
            }
        };

        static bool equal(const FT& a, const FT& b, FT epsilon = eps)
        {
            return a + epsilon > b && a < b + epsilon;
        }
        static bool less(const FT& a, const FT& b, FT epsilon = eps)
        {
            return a + epsilon <= b;
        }
        inline int ufs_find(std::vector<int>& f, size_t x) {
            if (f[x] == x) return x;
            return f[x] = ufs_find(f, f[x]);
        }
        void CalCenterLine(const CenterLineProps& centerlineProps);
        void connect_segments();
        void connect_segments2();

        Point_2 top_right;
        boost::shared_ptr<Ss> skeleton;
        std::vector<Location> locations;
        std::vector<PointData*> point_data_pool;
        std::vector<Segment_2> res_segments, sub_segments, smooth_segments;
        std::vector<std::pair<FT, FT>> res_seg_dis;
        std::vector<Point_2> log_points;
    };

    namespace SEH {
        using K = CGAL::Epick;
        using Ss = CGAL::Straight_skeleton_2<K>;
        using SsBuilderTraits = CGAL::Straight_skeleton_builder_traits_2<K>;
        using SsBuilder = CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss>;
        static void func(SsBuilder& ssb, boost::shared_ptr<Ss>& skeleton) {
            skeleton = ssb.construct_skeleton();
        }
        static void construct_skeleton(SsBuilder& ssb, boost::shared_ptr<Ss>& skeleton) {
            // Due to unknown CGAL issues, construct_skeleton() throws Exception Access Violation
            // This exception won't be caught by C++ standard Exception Handling (try...catch....)
            // To workaround this SEH exception, we wrap around the exception with MS extension __try...__except
            // When the SEH exception is caught, the std::exception is rethrew
            __try {
                func(ssb, skeleton);
            }
            __except (1) {
                throw("construct_skeleton exception ...\n");
            }
        }
    }
}

namespace CenterLine {
    template <typename K>
    Solver<K>::Solver(const std::string& geojson, const CenterLineProps& centerlineProps) {
        // set precision
        CGAL_setting();
        try {
            origin_space = GeoJSON::geojson_to_Pwh<K>(geojson);
            this->CalCenterLine(centerlineProps);
        }
        catch (const std::exception& e) {
            std::cerr << "!Error at SimplifyBoundary Module: " << e.what() << std::endl;
        }
        catch (...) {
            std::cerr << "!Caught unknown exception" << std::endl;
        }
    }

    template <typename K>
    void Solver<K>::CalCenterLine(const CenterLineProps& centerlineProps) {
        const Polygon_with_holes_2& polygon = origin_space;
        top_right = Point_2(polygon.outer_boundary().right_vertex()->x(), polygon.outer_boundary().top_vertex()->y());
        // total_edges is sum of the edges of outer boundary and holes.
        Polygon_2 outer = polygon.outer_boundary();
        size_t total_edges = outer.size(); // the size of line_data_pool.
        if (!polygon.outer_boundary().is_simple())
            throw  "Polygon is not simple!";

        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it) {
            if (!it->is_simple())
                throw  "Polygon is not simple!";
            total_edges += it->size();
        }


        // # get skeleton
        SsBuilder ssb;
        ssb.enter_contour(polygon.outer_boundary().vertices_begin(), polygon.outer_boundary().vertices_end());
        for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
            ssb.enter_contour(hole->vertices_begin(), hole->vertices_end());
        SEH::construct_skeleton(ssb, skeleton);

        // # unifon-find
        std::unordered_map<int, int> index; // vertex id -> index
        std::vector<int> f; // store the index of vertices according to the id, to use union-find set
        // std::cout << "size = " << skeleton->size_of_vertices() << std::endl;
        size_t cnt = 0;
        for (Vertex_iterator it = skeleton->vertices_begin(); it != skeleton->vertices_end(); ++it) {
            index[it->id()] = cnt;
            f.push_back(cnt);
            ++cnt;
        }
        // compress the skeleton vertices
        for (Halfedge_iterator it = skeleton->halfedges_begin(); it != skeleton->halfedges_end(); ++it)
            if (it->is_bisector()) {    // if the halfedge is not border edge, but skeleton edge
                Vertex_handle from = it->opposite()->vertex(), to = it->vertex();   // get two point of the edge
                if (equal(from->point().x(), to->point().x(), 1e-4) && equal(from->point().y(), to->point().y(), 1e-4)) {   // the location of a equals to b
                    int a = ufs_find(f, index[from->id()]), b = ufs_find(f, index[to->id()]);
                    f[a] = b;   // make a equal to b
                }
            }
        std::vector<int> rk(cnt, -1);
        int rk_it = 0;
        for (int i = 0; i < cnt; ++i) {
            f[i] = ufs_find(f, i);
            if (rk[f[i]] == -1) rk[f[i]] = rk_it++;
            f[i] = rk[f[i]];
        }

        // # store information of point and time in locations
        std::unordered_set<std::pair<int, int>, PairHash> edge_set; // store edges that have already been processed.
        locations.resize(rk_it);    // rk_it is number of connected components
        for (Vertex_iterator it = skeleton->vertices_begin(); it != skeleton->vertices_end(); ++it) {
            int id = f[index[it->id()]];
            locations[id].point = it->point();
            //locations[cnt].time = it->time();
            locations[id].time = it->time() * it->time();   // time of contour vertex is 0
        }

        FT MIN_SQUARED_COS_THETA = std::cos(centerlineProps.theta_thre3 * M_PI / 180);
        MIN_SQUARED_COS_THETA *= MIN_SQUARED_COS_THETA;
        FT THINK_CLOSE_PERP_MIN = std::cos(120 * M_PI / 180);
        FT THINK_CLOSE_PERP_MAX = std::cos(60 * M_PI / 180);
        for (Halfedge_iterator it = skeleton->halfedges_begin(); it != skeleton->halfedges_end(); ++it)
            if (it->is_bisector()) {    // the halfedge is not border
                Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
                double len_ = CGAL::squared_distance(from->point(), to->point());
                int a = f[index[from->id()]], b = f[index[to->id()]];
                Halfedge_handle l = it->defining_contour_edge(), r = it->opposite()->defining_contour_edge();   // get the two contour edges of the skeleton edge
                /*
                      Halfedge_handle  defining_contour_edge() { return this->face()->halfedge() ; }
                      this->face()->halfedge(): return the first halfedge of this face, which is the halfedge of the contour, because others are generated by it.
                */
                if (a < b) {
                    if (from->time() > to->time()) {
                        std::swap(a, b);
                        std::swap(from, to);
                        std::swap(l, r);
                    }
                    if (edge_set.count(std::make_pair(a, b))) continue;
                    edge_set.insert(std::make_pair(a, b));

                    // # save skeleton edge to edge_data_pool, inlude all skeleton edge
                    int e_id = point_data_pool.size();
                    //PointData *edge = new PointData(locations, index[from->id()], index[to->id()], e_id);
                    PointData* edge = new PointData(locations, a, b, e_id);
                    point_data_pool.push_back(edge);

                    // # filter the edges

                    //std::cout << "edge = " << from->point() << " " << to->point() << std::endl;
                    Point_2 v_la = l->vertex()->point();
                    Point_2 v_lb = l->opposite()->vertex()->point();
                    Vector_2 v_l = v_la - v_lb;
                    Point_2 v_ra = r->vertex()->point();
                    Point_2 v_rb = r->opposite()->vertex()->point();
                    Vector_2 v_r = v_ra - v_rb;
                    // std::cout << "l = " << l->opposite()->vertex()->point() << " " << l->vertex()->point() << std::endl;
                    // std::cout << "r = " << r->opposite()->vertex()->point() << " " << r->vertex()->point() << std::endl;
                    // std::cout << "v_l = " << v_l << "\nv_r = " << v_r << std::endl;
                    FT inner = v_l.x() * v_r.x() + v_l.y() * v_r.y();
                    FT outer = v_l.x() * v_r.y() - v_l.y() * v_r.x();
                    FT absolute_scale = v_l.squared_length() * v_r.squared_length();
                    FT squared_cos_theta = (inner * inner) / absolute_scale;
                    edge->squared_cos_theta = squared_cos_theta;
                    edge->inner = inner;
                    //edge->is_ans = (inner < 0 && outer >= 0 && inner * inner * 4 >= absolute_scale);
                    //edge->is_ans = (!equal(v_la.x(), v_rb.x()) || !equal(v_la.y(), v_rb.y()));
                    edge->is_ans = !(equal(v_la.x(), v_rb.x()) && equal(v_la.y(), v_rb.y()) && squared_cos_theta >= MIN_SQUARED_COS_THETA && inner >= 0);
                    //edge->is_ans = (inner < 0 && outer >= 0 && squared_cos_theta >= MIN_SQUARED_COS_THETA);
                    //edge->is_ans = inner < 0;
                    //edge->is_ans = (inner < 0 && outer >= 0);
                    // edge->is_ans = (outer >= 0);
                    /*
                        inner < 0: The angel of two border vectors is greater than 90 degrees.
                        outer >= 0: The two border vectors are counter-clockwise direction.
                        inner * inner * 4 >= absolute_scale: The angel of two border vectors not in (90-theta, 90+theta)
                    */

                    // std::cout << edge->is_ans << std::endl;
                    //locations[index[from->id()]].branches.emplace_back(edge, 0);
                    //locations[index[to->id()]].branches.emplace_back(edge, 1);
                    locations[a].branches.emplace_back(edge, 0);
                    locations[b].branches.emplace_back(edge, 1);
                }
            }

        // find the leaf branch
        FT outer_width = outer.right_vertex()->x() - outer.left_vertex()->x();
        FT outer_height = outer.top_vertex()->y() - outer.bottom_vertex()->y();
        FT MIN_CONVEX_HEIGHT = CGAL::sqrt(outer_width * outer_width + outer_height * outer_height) * centerlineProps.min_convex_height;
        std::unordered_set<int> stop_loc_ids;
        std::unordered_set< PointData*> start_edges;
        for (PointData* e : point_data_pool)
            if (e->is_ans) {
                int ans_num1 = locations[e->start_loc].get_ans_num();
                if (ans_num1 == 1)
                    start_edges.insert(e);
                int ans_num2 = locations[e->end_loc].get_ans_num();
                if (e->inner < 0 && e->squared_cos_theta >= MIN_SQUARED_COS_THETA) {
                    stop_loc_ids.insert(e->end_loc);
                    stop_loc_ids.insert(e->start_loc);
                }
            }
        for(auto e : start_edges)
            log_points.push_back(locations[e->start_loc].point);
        for (PointData* e : start_edges) {
            PointData* nowEdge = e;
            while (nowEdge != nullptr && stop_loc_ids.find(nowEdge->start_loc) == stop_loc_ids.end()) {
                nowEdge->is_ans = false;
                nowEdge = locations[nowEdge->end_loc].get_next_ans_edge();
            }

            /*while (nowEdge != nullptr && stop_loc_ids.find(nowEdge->start_loc) == stop_loc_ids.end()) {
                nowEdge->is_ans = false;
                if(locations[e->start_loc].get_ans_num() == 1)
                    nowEdge = locations[nowEdge->start_loc].get_next_ans_edge();
                else
                    nowEdge = locations[nowEdge->end_loc].get_next_ans_edge();
            }

            if (nowEdge != nullptr) {
                FT lenEdge = CGAL::sqrt(CGAL::squared_distance(nowEdge->source(), nowEdge->target()));
                nowEdge->is_ans = false;
                PointData* nextEdge;
                if (locations[nowEdge->start_loc].get_ans_num() == 1)
                    nextEdge = locations[nowEdge->start_loc].get_next_ans_edge();
                else
                    nextEdge = locations[nowEdge->end_loc].get_next_ans_edge();
                if (nextEdge == nullptr) {
                    nowEdge->is_ans = true;
                    continue;
                }
                if (lenEdge <= MIN_CONVEX_HEIGHT && nextEdge->squared_cos_theta >= THINK_CLOSE_PERP_MIN && nextEdge->squared_cos_theta <= THINK_CLOSE_PERP_MAX) {
                    nextEdge->is_ans = false;
                    nowEdge = locations[nextEdge->end_loc].get_next_ans_edge();
                }
                else {
                    nowEdge->is_ans = true;
                    continue;
                }

                while (nowEdge != nullptr && stop_loc_ids.find(nowEdge->start_loc) == stop_loc_ids.end()) {
                    nowEdge->is_ans = false;
                    nowEdge = locations[nowEdge->end_loc].get_next_ans_edge();
                }
            }*/
        }

        // # connect segments
        size_t sep = point_data_pool.size();
        try {
            this->connect_segments2();
        }
        catch (const char* str) {
            std::cerr << "connect segments error: " << str << std::endl;
        }

        if (centerlineProps.isSmooth == false) {
            // # get res_segments, res_seg_dis and sub_segments
            int counter = 0;
            for (PointData* it : point_data_pool) {
                ++counter;
                if (counter - 1 == sep) this->seg_cnt_before_connect = this->res_segments.size();
                if (it->end_loc == -1) {
                    std::cout << "point " << counter - 1 << "not ended" << std::endl;
                    continue;
                }
                if (it->is_ans) {
                    if (it->end_loc != -1) {
                        this->res_segments.emplace_back(it->source(), it->target());
                        this->res_seg_dis.emplace_back(it->start_time(), it->end_time());
                    }
                    else {
                        std::cerr << "point " << it->point_id << " not ended" << std::endl;
                    }
                }
                else
                    this->sub_segments.emplace_back(it->source(), it->target());
            }
            if (point_data_pool.size() == sep) this->seg_cnt_before_connect = this->res_segments.size();
        }
        else {
            std::unordered_map<int, int> loc2degree;
            int counter = 0;
            for (PointData* it : point_data_pool) {
                ++counter;
                if (it->end_loc == -1) {
                    std::cout << "point " << counter - 1 << "not ended" << std::endl;
                    continue;
                }
                if (it->is_ans) {
                    if (it->end_loc != -1) {
                        if (loc2degree.find(it->start_loc) == loc2degree.end())
                            loc2degree[it->start_loc] = 1;
                        else
                            loc2degree[it->start_loc]++;
                        if (loc2degree.find(it->end_loc) == loc2degree.end())
                            loc2degree[it->end_loc] = 1;
                        else
                            loc2degree[it->end_loc]++;
                    }
                    else {
                        std::cerr << "point " << it->point_id << " not ended" << std::endl;
                    }
                }
            }

            // TODO: fix out smooth_segment2
            //smooth_segments2(loc2degree, smooth_segments, Sprops, centerlineProps.insertNum);
        }
    }

    template <typename K>
    inline void Solver<K>::connect_segments() {
        for (int i = 0; i < locations.size(); ++i) {
            std::cout << "location " << i << ": <" << locations[i].point << ", " << locations[i].time << ">" << std::endl;
        }
        for (auto p : point_data_pool) {
            std::cout << p->point_id << ":" << std::endl;
            std::cout << "start_loc = " << p->start_loc << std::endl;
            std::cout << "end_loc = " << p->end_loc << std::endl;
        }

        std::vector<bool> visited(locations.size(), 0), inQ(locations.size(), 0);
        std::vector<std::vector<Vector_2>> in_vectors(locations.size());
        std::priority_queue<std::pair<FT, int>> Q; // <-time, loc_id>

        for (auto p : point_data_pool) if (p->is_ans) {
            if (p->end_loc == -1) {
                std::cerr << "point_data " << p->point_id << " not ended" << std::endl;
            }
            else {
                in_vectors[p->start_loc].push_back(p->point(0) - p->point(1));
                in_vectors[p->end_loc].push_back(p->point(1) - p->point(0));
            }
        }
        std::cout << "preprocessing finished" << std::endl;
        size_t locations_size = locations.size();
        for (int i = 0; i < locations.size(); ++i) if (locations[i].branches.size() == 1) { // leaves
            visited[i] = true;
            auto& e = locations[i].branches[0];
            if (e.first->end_loc == -1) continue;
            int new_loc = e.first->location(e.second ^ 1);
            std::cout << "i=" << i << "\te.first->point_id=" << e.first->point_id << std::endl;
            std::cout << "e.second = " << e.second << std::endl;
            std::cout << "new_loc = " << new_loc << std::endl;
            if (e.first->is_ans) {
                Vector_2 v = e.first->target() - e.first->source();
                if (e.second == 1) v = -v;
                std::cout << v;
                in_vectors[new_loc].push_back(v);
            }
            if (!inQ[new_loc]) {
                inQ[new_loc] = true;
                Q.emplace(-locations[new_loc].time, new_loc);
            }
        }
        while (!Q.empty()) {
            int loc_id = Q.top().second; Q.pop();
            visited[loc_id] = true;

            for (size_t i = 0; i < locations[loc_id].branches.size(); ++i) {
                const auto e = locations[loc_id].branches[i];
                if (e.second == false) {
                    if (e.first->end_loc == -1) {
                        std::cerr << "Point " << e.first->point_id << " not ended." << std::endl;
                        continue;
                    }
                    int new_loc = e.first->location(e.second ^ 1);
                    std::cout << "\ncur = " << loc_id << "\tnew_loc=" << new_loc << std::endl;
                    if (visited[new_loc]) continue;
                    if (!e.first->is_ans && in_vectors[loc_id].size() != 0) {
                        const Point_2& src = e.first->point(e.second), & dest = e.first->point(e.second ^ 1);
                        Vector_2 base_v = dest - src;
                        std::cout << "connecting PointData " << e.first->point_id << std::endl;
                        std::cout << "from: "; for (auto v : in_vectors[loc_id]) std::cout << v << " "; std::cout << std::endl;
                        std::cout << e.first->point(e.second) << " => " << e.first->point(e.second ^ 1) << std::endl;
                        std::cout << "vect = " << e.first->point(e.second ^ 1) - e.first->point(e.second) << std::endl;
                        std::cout << "to: "; for (auto v : in_vectors[new_loc]) std::cout << -v << " "; std::cout << std::endl;

                        if (equal(src.x(), dest.x()) && equal(src.y(), dest.y())) { // �˻��ĵ�
                            throw("repeated point");
                            //in_vectors[new_loc].insert(in_vectors[new_loc].end(), in_vectors[loc_id].begin(), in_vectors[loc_id].end());
                        }
                        else if (!_context.prefer_ortho) {
                            PointData* new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                            point_data_pool.push_back(new_point);
                            in_vectors[loc_id].push_back(-base_v);
                            in_vectors[new_loc].push_back(base_v);
                        }
                        else if (in_vectors[new_loc].size() == 0) { // ��ص�����PointData�����Ǵ�
                            Vector_2 used_vector;
                            FT value = -1;
                            for (auto v : in_vectors[loc_id]) {
                                evalDirection(v, base_v, used_vector, value);
                            }
                            if (value == 1 || value < 0) {
                                PointData* new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                                point_data_pool.push_back(new_point);
                                in_vectors[loc_id].push_back(-base_v);
                                in_vectors[new_loc].push_back(base_v);
                            }
                            else {
                                Point_2 foot = Line_2(src, used_vector).projection(dest);
                                int turning = locations.size();
                                locations.push_back(Location(foot, (locations[loc_id].time + locations[new_loc].time) / 2));
                                PointData* point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                                point_data_pool.push_back(point0);
                                in_vectors[loc_id].push_back(src - foot);

                                PointData* point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                                point_data_pool.push_back(point1);
                                in_vectors[new_loc].push_back(dest - foot);
                            }
                        }
                        else { // �������in_vectors
                            Vector_2 used_vector0, used_vector1, used_v0, used_v1;
                            FT value = -2, cur_value;
                            //bool is_perp, use_perp = false;
                            //Vector_2 perp_v0, perp_v1;
                            for (auto v0 : in_vectors[loc_id]) {
                                FT Sin0 = outer_product(v0, base_v) / CGAL::sqrt(v0.squared_length() * base_v.squared_length());
                                FT Cos0 = Cosine(v0, base_v);
                                for (auto v1 : in_vectors[new_loc]) {
                                    FT Sin1 = outer_product(base_v, -v1) / CGAL::sqrt(v1.squared_length() * base_v.squared_length());
                                    FT Cos1 = Cosine(-v1, base_v);
                                    //is_perp = false;
                                    if (equal(Cos0, 1) || equal(Cos1, 1)) {
                                        used_v0 = used_v1 = base_v;
                                        cur_value = 4;
                                    }
                                    //else if(Cos0 < 0 || Cos1 < 0){
                                    //    used_v0 = used_v1 = base_v;
                                    //    if(Cos0 >= 0) Cos0 = 1 + (Cos0 - 1) * Cos0;
                                    //    if(Cos1 >= 0) Cos1 = 1 + (Cos1 - 1) * Cos1;
                                    //    cur_value = (Cos0 + Cos1) / 2;
                                    //}
                                    else if (Cos0 <= 0 || Cos1 <= 0) {
                                        Vector_2 nv0 = v0, nv1 = v1;
                                        //is_perp = true;
                                        if (Cos0 <= 0) {
                                            FT tmp;
                                            if (Sin0 < 0) {
                                                nv0 = v0.perpendicular(CGAL::CLOCKWISE);
                                                tmp = -Sin0; Sin0 = Cos0;
                                            }
                                            else {
                                                nv0 = v0.perpendicular(CGAL::COUNTERCLOCKWISE);
                                                tmp = Sin0; Sin0 = -Cos0;
                                            }
                                            Cos0 = tmp;
                                        }
                                        if (Cos1 <= 0) {
                                            FT tmp;
                                            if (Sin1 < 0) {
                                                nv1 = v1.perpendicular(CGAL::COUNTERCLOCKWISE);
                                                tmp = -Sin1; Sin1 = Cos1;
                                            }
                                            else {
                                                nv1 = v1.perpendicular(CGAL::CLOCKWISE);
                                                tmp = Sin1; Sin1 = -Cos1;
                                            }
                                            Cos1 = tmp;
                                        }
                                        if (equal(Cos0, 1) || equal(Cos1, 1)) {
                                            used_v0 = used_v1 = base_v;
                                            cur_value = 2;
                                        }
                                        else if (Sin0 * Sin1 < 0) {
                                            used_v0 = nv0; used_v1 = -nv1;
                                            FT tmp = Cosine(nv0, -nv1);
                                            cur_value = tmp * tmp / 2;
                                        }
                                        else {
                                            //used_v0 = nv0; used_v1 = -nv1;
                                            //FT t = Sine(nv0, -nv1);
                                            //cur_value = t * t;
                                            FT inner = nv0.x() * nv1.x() + nv0.y() * nv1.y();
                                            FT absol = nv0.squared_length() * nv1.squared_length();
                                            if (inner > 0 && inner * inner * 4 >= absol) { // >= 120 degrees
                                                used_v0 = used_v1 = base_v;
                                                FT t0 = Cos0 * Cos0 * 2 - 1; t0 = t0 * t0;
                                                FT t1 = Cos1 * Cos1 * 2 - 1; t1 = t1 * t1;
                                                cur_value = (t0 + t1) / 4;
                                            }
                                            else {
                                                used_v0 = nv0; used_v1 = -nv1;
                                                FT tmp = Sine(nv0, -nv1);
                                                cur_value = tmp * tmp;
                                            }
                                        }
                                    }
                                    else if (Sin0 * Sin1 < 0) {
                                        used_v0 = v0; used_v1 = -v1;
                                        //cur_value = -Sin0 * Sin1;
                                        FT tmp = Cosine(v0, -v1);
                                        cur_value = tmp * tmp;
                                    }
                                    else { // Sin0 * Sin1 > 0
                                        FT inner = v0.x() * v1.x() + v0.y() * v1.y();
                                        FT absol = v0.squared_length() * v1.squared_length();
                                        if (inner > 0 && inner * inner * 4 >= absol) { // >= 120 degrees
                                            used_v0 = used_v1 = base_v;
                                            FT t0 = Cos0 * Cos0 * 2 - 1; t0 = t0 * t0;
                                            FT t1 = Cos1 * Cos1 * 2 - 1; t1 = t1 * t1;
                                            cur_value = (t0 + t1) / 2;
                                        }
                                        else {
                                            used_v0 = v0; used_v1 = -v1;
                                            FT tmp = Sine(v0, -v1);
                                            cur_value = 2 * tmp * tmp;
                                        }
                                    }
                                    if (cur_value > value) {
                                        //if(is_perp){ use_perp = true; perp_v0 = v0; perp_v1 = v1; }
                                        //else use_perp = false;
                                        value = cur_value;
                                        used_vector0 = used_v0, used_vector1 = used_v1;
                                    }
                                }
                            }
                            //if(use_perp){
                            //    Vector_2 v0 = perp_v0, v1 = perp_v1; Vector_2 nv0 = v0, nv1 = v1;
                            //    FT Sin0 = outer_product(v0, base_v) / CGAL::sqrt(v0.squared_length() * base_v.squared_length());
                            //    FT Sin1 = outer_product(base_v, -v1) / CGAL::sqrt(v1.squared_length() * base_v.squared_length());
                            //    FT Cos0 = Cosine(v0, base_v); FT Cos1 = Cosine(-v1, base_v);
                            //    std::cout << "use_perp: " << src << " "  << dest << std::endl;
                            //    std::cout << v0 << " " << base_v << " " << v1 << std::endl;
                            //    //...
                            //}
                            // end DEBUG
                            if (used_vector0 == base_v && used_vector1 == base_v) { // ֱ������
                                PointData* new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                                point_data_pool.push_back(new_point);
                                in_vectors[loc_id].push_back(-base_v);
                                in_vectors[new_loc].push_back(base_v);
                            }
                            else if (outer_product(used_vector0, base_v) * outer_product(base_v, used_vector1) < 0) {
                                Vector_2 v0 = used_vector0 / CGAL::sqrt(used_vector0.squared_length());
                                Vector_2 v1 = used_vector1 / CGAL::sqrt(used_vector1.squared_length());
                                Point_2 mid = src + (dest - src) / 2, foot0, foot1;
                                if (inner_product(v0, base_v) < inner_product(v1, base_v)) {
                                    foot0 = Line_2(src, used_vector0).projection(mid);
                                    Line_2 mid_line(mid, used_vector0.perpendicular(CGAL::LEFT_TURN));
                                    auto inter = CGAL::intersection(mid_line, Line_2(dest, used_vector1));
                                    Point_2* tmp_p;
                                    if (!inter || !(tmp_p = boost::get<Point_2>(&*inter))) throw("no intersection of mid_line and v1");
                                    foot1 = *tmp_p;
                                }
                                else {
                                    foot1 = Line_2(dest, used_vector1).projection(mid);
                                    Line_2 mid_line(mid, used_vector1.perpendicular(CGAL::LEFT_TURN));
                                    auto inter = CGAL::intersection(mid_line, Line_2(src, used_vector0));
                                    Point_2* tmp_p;
                                    if (!inter || !(tmp_p = boost::get<Point_2>(&*inter))) throw("no intersection of mid_line and v0");
                                    foot0 = *tmp_p;
                                }
                                int turning0 = locations.size();
                                locations.push_back(Location(foot0, (locations[loc_id].time * 3 + locations[new_loc].time) / 4));
                                int turning1 = locations.size();
                                locations.push_back(Location(foot1, (locations[loc_id].time + locations[new_loc].time * 3) / 4));
                                PointData* point0 = new PointData(locations, loc_id, turning0, point_data_pool.size());
                                point_data_pool.push_back(point0);
                                in_vectors[loc_id].push_back(src - foot0);

                                PointData* point1 = new PointData(locations, turning0, turning1, point_data_pool.size());
                                point_data_pool.push_back(point1);
                                PointData* point2 = new PointData(locations, turning1, new_loc, point_data_pool.size());
                                point_data_pool.push_back(point2);
                                in_vectors[new_loc].push_back(dest - foot1);
                            }
                            else {
                                auto inter = CGAL::intersection(Line_2(src, used_vector0), Line_2(dest, used_vector1));
                                Point_2* tmp_p;
                                if (!inter || !(tmp_p = boost::get<Point_2>(&*inter))) throw("no intersection of v0 and v1");
                                int turning = locations.size();
                                locations.push_back(Location(*tmp_p, (locations[loc_id].time + locations[new_loc].time) / 2));
                                PointData* point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                                point_data_pool.push_back(point0);
                                in_vectors[loc_id].push_back(src - *tmp_p);

                                PointData* point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                                point_data_pool.push_back(point1);
                                in_vectors[new_loc].push_back(dest - *tmp_p);
                            }
                        }
                    }
                    else if (in_vectors[loc_id].size() == 0 && locations[loc_id].branches.size() > 3) {
                        const Point_2& src = e.first->point(e.second), & dest = e.first->point(e.second ^ 1);
                        Vector_2 base_v = dest - src;
                        if (in_vectors[new_loc].size()) {
                            Vector_2 used_vector;
                            FT value = -1;
                            for (auto v : in_vectors[new_loc]) {
                                evalDirection(v, -base_v, used_vector, value);
                            }
                            if (value == 1 || value < 0) { // simply link
                                PointData* new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                                point_data_pool.push_back(new_point);
                                in_vectors[new_loc].push_back(base_v);
                                in_vectors[loc_id].push_back(-base_v);
                            }
                            else {
                                Point_2 foot = Line_2(src, used_vector).projection(dest);
                                int turning = locations.size();
                                locations.push_back(Location(foot, (locations[loc_id].time + locations[new_loc].time) / 2));
                                PointData* point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                                point_data_pool.push_back(point0);
                                in_vectors[loc_id].push_back(src - foot);

                                PointData* point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                                point_data_pool.push_back(point1);
                                in_vectors[new_loc].push_back(dest - foot);
                            }
                        }
                        else { // in_vectors of src and dest are both empty
                            if (equal(base_v.x(), 0) || equal(base_v.y(), 0)) {
                                PointData* new_point = new PointData(locations, loc_id, new_loc, point_data_pool.size());
                                point_data_pool.push_back(new_point);
                                in_vectors[new_loc].push_back(base_v);
                                in_vectors[loc_id].push_back(-base_v);
                            }
                            else {
                                Point_2 foot;
                                if (base_v.x() < base_v.y()) foot = Point_2(src.x(), dest.y());
                                else foot = Point_2(dest.x(), src.y());
                                int turning = locations.size();
                                locations.push_back(Location(foot, (locations[loc_id].time + locations[new_loc].time) / 2));
                                PointData* point0 = new PointData(locations, loc_id, turning, point_data_pool.size());
                                point_data_pool.push_back(point0);
                                in_vectors[loc_id].push_back(src - foot);

                                PointData* point1 = new PointData(locations, turning, new_loc, point_data_pool.size());
                                point_data_pool.push_back(point1);
                                in_vectors[new_loc].push_back(dest - foot);
                            }
                        }
                    }
                    if (!inQ[new_loc]) {
                        inQ[new_loc] = true;
                        Q.emplace(-locations[new_loc].time, new_loc);
                    }
                }
            }
        }
        // extend leaf segments
        std::vector<int> degree(locations.size(), 0);
        for (auto p : point_data_pool) if (p->is_ans) {
            ++degree[p->start_loc];
            std::cout << p->point_id << " is ans\n";
            std::cout << "ans=" << p->source() << " " << p->target() << std::endl;
            if (p->end_loc != -1) ++degree[p->end_loc];
        }
        for (size_t loc_id = 0; loc_id != locations_size; ++loc_id) {
            //if(degree[loc_id] == 1){
            if (degree[loc_id] == 1) {
                if (in_vectors[loc_id].size() == 0) {
                    throw("connect_segments error: leaf node has no in_vector");
                }
                std::cout << "degree" << in_vectors[loc_id].size() << std::endl;
                std::cout << "leaf point " << loc_id << ": " << locations[loc_id].point << std::endl;
                std::cout << "vector=" << in_vectors[loc_id][0] << std::endl;
                try {
                    Point_2 new_p = get_intersection(*origin_space, locations[loc_id].point, in_vectors[loc_id][0]);
                    if (CGAL::squared_distance(new_p, locations[loc_id].point) < 1e-6) {
                        continue;
                    }
                    std::cout << "intersection " << locations[loc_id].point << " " << new_p << std::endl;
                    size_t endpoint = locations.size();
                    locations.push_back(Location(new_p, locations[loc_id].time));
                    PointData* point0 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                    point_data_pool.push_back(point0);
                }
                catch (const char* str) {
                    std::cerr << "connect_segments error = \n" << str << std::endl;
                }
            }
            else if (degree[loc_id] == 0 && locations[loc_id].branches.size() > 2) {
                bool has_out = false;
                for (auto& branch : locations[loc_id].branches) if (branch.second == 0) {
                    has_out = true;
                    break;
                }
                if (!has_out) {
                    std::cout << "isoloated point" << locations[loc_id].point << std::endl;
                    try {
                        Point_2 l = get_intersection(*origin_space, locations[loc_id].point, Vector_2(-1, 0));
                        std::cout << "intersection " << locations[loc_id].point << " " << l << std::endl;
                        size_t endpoint = locations.size();
                        locations.push_back(Location(l, locations[loc_id].time));
                        PointData* point0 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                        point_data_pool.push_back(point0);

                        Point_2 r = get_intersection(*origin_space, locations[loc_id].point, Vector_2(1, 0));
                        std::cout << "intersection " << locations[loc_id].point << " " << r << std::endl;
                        endpoint = locations.size();
                        locations.push_back(Location(r, locations[loc_id].time));
                        PointData* point1 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                        point_data_pool.push_back(point1);
                    }
                    catch (const char* str) {
                        std::cerr << "connect_segments error = \n" << str << std::endl;
                    }
                }
            }
        }
    }

    template <typename K>
    inline void Solver<K>::connect_segments2() {
        // # init
        std::vector<bool> visited(locations.size(), 0), inQ(locations.size(), 0);
        std::vector<std::vector<Vector_2>> in_vectors(locations.size());    // the in vectors which are the true answer
        std::priority_queue<std::pair<FT, int>> Q; // <-time, loc_id>

        // # preprocess
        for (auto p : point_data_pool)
            if (p->is_ans) {
                if (p->end_loc == -1) {
                    std::cerr << "point_data " << p->point_id << " not ended" << std::endl;
                }
                else {
                    in_vectors[p->start_loc].push_back(p->point(0) - p->point(1));  // locations[start_loc].point - locations[end_loc].point
                    in_vectors[p->end_loc].push_back(p->point(1) - p->point(0));    // locations[end_loc].point - locations[start_loc].point
                }
            }
        std::cout << "preprocessing finished" << std::endl;

        // # leaf point
        size_t locations_size = locations.size();
        for (int i = 0; i < locations.size(); ++i)
            if (locations[i].branches.size() == 1) { // leaves
                visited[i] = true;
                auto& e = locations[i].branches[0];
                if (e.first->end_loc == -1) continue;
                int new_loc = e.first->location(e.second ^ 1);
                /*std::cout << "i=" << i << "\te.first->point_id=" << e.first->point_id << std::endl;
                std::cout << "e.second = " << e.second << std::endl;
                std::cout << "new_loc = " << new_loc << std::endl;*/
                if (e.first->is_ans) {
                    Vector_2 v = e.first->target() - e.first->source();
                    if (e.second == 1) v = -v;
                    // std::cout << v;
                    in_vectors[new_loc].push_back(v);
                }
                if (!inQ[new_loc]) {
                    inQ[new_loc] = true;
                    Q.emplace(-locations[new_loc].time, new_loc);
                }
            }

        // # extend leaf segments
        std::vector<int> degree(locations.size(), 0);
        for (auto p : point_data_pool) if (p->is_ans) {
            ++degree[p->start_loc];
            // std::cout << p->point_id << " is ans\n";
            // std::cout << "ans=" << p->source() << " " << p->target() << std::endl;
            if (p->end_loc != -1) ++degree[p->end_loc];
        }
        for (size_t loc_id = 0; loc_id != locations_size; ++loc_id) {
            //if(degree[loc_id] == 1){
            if (degree[loc_id] == 1) {
                if (in_vectors[loc_id].size() == 0) {
                    throw("connect_segments error: leaf node has no in_vector");
                }
                //std::cout << "degree" << in_vectors[loc_id].size() << std::endl;
                //std::cout << "leaf point " << loc_id << ": " << locations[loc_id].point << std::endl;
                //std::cout << "vector=" << in_vectors[loc_id][0] << std::endl;
                try {
                    Point_2 new_p = get_intersection2<K>(origin_space, locations[loc_id].point, in_vectors[loc_id][0]);
                    if (CGAL::squared_distance(new_p, locations[loc_id].point) < 1e-6) {
                        continue;
                    }
                    // std::cout << "intersection " << locations[loc_id].point << " " << new_p << std::endl;
                    size_t endpoint = locations.size();
                    locations.push_back(Location(new_p, locations[loc_id].time));
                    PointData* point0 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                    point_data_pool.push_back(point0);
                    locations[endpoint].branches.push_back(std::make_pair(point0, 1));
                    locations[loc_id].branches.push_back(std::make_pair(point0, 0));
                }
                catch (const char* str) {
                    std::cerr << "connect_segments error = \n" << str << std::endl;
                }
            }
            else if (degree[loc_id] == 0 && locations[loc_id].branches.size() > 2) {
                bool has_out = false;
                for (auto& branch : locations[loc_id].branches) if (branch.second == 0) {
                    has_out = true;
                    break;
                }
                if (!has_out) {
                    //std::cout << "isoloated point" << locations[loc_id].point << std::endl;
                    try {
                        Point_2 l = get_intersection2<K>(origin_space, locations[loc_id].point, Vector_2(-1, 0));
                        //std::cout << "intersection " << locations[loc_id].point << " " << l << std::endl;
                        size_t endpoint = locations.size();
                        locations.push_back(Location(l, locations[loc_id].time));
                        PointData* point0 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                        point_data_pool.push_back(point0);

                        Point_2 r = get_intersection2<K>(origin_space, locations[loc_id].point, Vector_2(1, 0));
                        //std::cout << "intersection " << locations[loc_id].point << " " << r << std::endl;
                        endpoint = locations.size();
                        locations.push_back(Location(r, locations[loc_id].time));
                        PointData* point1 = new PointData(locations, loc_id, endpoint, point_data_pool.size());
                        point_data_pool.push_back(point1);
                    }
                    catch (const char* str) {
                        std::cerr << "connect_segments error = \n" << str << std::endl;
                    }
                }
            }
        }
    }
}
#endif // CENTERLINE_HPP