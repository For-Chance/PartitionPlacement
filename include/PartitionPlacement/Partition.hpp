#ifndef PARTITION_HPP
#define PARTITION_HPP
#include "stdafx.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <list>
#include <chrono>

#include "SimplifyBoundary.hpp"

namespace Partition
{
    template <typename K>
    struct PartitionProps {
        bool withSimplifyBoundary;
        PartitionProps() {
            withSimplifyBoundary = true;
        }
    };

    template <typename K>
    class Solver
    {
        using FT = typename K::FT;
        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;
		using Segment_2 = CGAL::Segment_2<K>;
        using Polygon_list = std::list<Polygon_2>;
        using Vector_2 = CGAL::Vector_2<K>;
        using Line_2 = CGAL::Line_2<K>;
        using Ss = CGAL::Straight_skeleton_2<K>;
        using Halfedge_iterator = typename Ss::Halfedge_iterator;
        using Halfedge_handle = typename Ss::Halfedge_handle;
        using Vertex_handle = typename Ss::Vertex_handle;
        using Vertex_iterator = typename Ss::Vertex_iterator;
		using Face_handle = typename Ss::Face_handle;
		using Face_iterator = typename Ss::Face_iterator;
        using SsBuilderTraits = CGAL::Straight_skeleton_builder_traits_2<K>;
        using SsBuilder = CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss>;

        using PartitionProps = PartitionProps<K>;
    private:
        /// <summary>
		/// build skeleton of a Polygon_with_holes_2
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        boost::shared_ptr<Ss> build_skeleton(const Polygon_with_holes_2& polygon) {
            std::cout << "start building skeleton..." << std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            SsBuilder ssb;
            ssb.enter_contour(polygon.outer_boundary().vertices_begin(), polygon.outer_boundary().vertices_end());
            for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole)
                ssb.enter_contour(hole->vertices_begin(), hole->vertices_end());
			boost::shared_ptr<Ss> skeleton = ssb.construct_skeleton();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
			std::cout << "skeleton built! Time elapsed:" << diff.count() << " s, num of skeleton face:" << skeleton->size_of_faces()<< std::endl;
            return skeleton;
        };
        /// <summary>
        /// get  Polygon_2 from a Face_handle
        /// </summary>
        /// <param name="face"></param>
        /// <returns></returns>
        Polygon_2 get_Poly_from_Face(const Face_handle& face) {
            Polygon_2 poly;
            Halfedge_handle h = face->halfedge();
            Halfedge_handle first = h;
            do {
                Point_2 point = h->vertex()->point();
                poly.push_back(point);
                h = h->next();
            } while (h != first);
            return poly;
        }
    public:
        Solver() {}
        Solver(const Polygon_with_holes_2& space, const PartitionProps& props = PartitionProps());
        
        Polygon_with_holes_2 origin_space;
        Polygon_with_holes_2 polygon;
        boost::shared_ptr<Ss> skeleton;
		std::vector<Segment_2> skeleton_segments;
        std::vector<Segment_2> skeleton_centerlines;
        std::vector<Segment_2> skeleton_otherlines;
        std::vector<Polygon_2> skeleton_faces;
        std::vector<std::vector<Polygon_2>> init_partition;
        std::vector<Polygon_2> uncertain_parts;
        std::vector<Segment_2> split_segments;
    };
}

namespace Partition
{
    template <typename K>
    Solver<K>::Solver(const Polygon_with_holes_2& space, const PartitionProps& props)
    {
		this->origin_space = space;
        if (props.withSimplifyBoundary) {
            // 1. Simplify Boundary
			SimplifyBoundary::Solver<K> sbSolver(this->origin_space);
            this->polygon = sbSolver.simplify_space;
            
            // 2. skeleton 
			this->skeleton = build_skeleton(this->polygon);
            for (Halfedge_iterator it = this->skeleton->halfedges_begin(); it != this->skeleton->halfedges_end(); ++it) {
                Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
				skeleton_segments.push_back(Segment_2(from->point(), to->point()));
            }
            
			std::cout << "partition..." << std::endl;
            // 3. partition
            // 3.1 save to skeleton_faces
            for (Face_iterator it = this->skeleton->faces_begin(); it != this->skeleton->faces_end(); ++it) {
                Face_handle face = it;
                skeleton_faces.push_back(get_Poly_from_Face(face));
            }
            // 3.2 find center line
            // first we try to make skeleton of no connect boundary is center line
            std::unordered_map<Vertex_handle, int> centerline_vertex2cnt;
            for (Halfedge_iterator it = this->skeleton->halfedges_begin(); it != this->skeleton->halfedges_end(); ++it) {
                if (it->is_bisector()) {    // the halfedge is not border
                    Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
                    Halfedge_handle l = it->defining_contour_edge(), r = it->opposite()->defining_contour_edge();   // get the two contour edges of the skeleton edge
                    /*
                          Halfedge_handle  defining_contour_edge() { return this->face()->halfedge() ; }
                          this->face()->halfedge(): return the first halfedge of this face, which is the halfedge of the contour, because others are generated by it.
                    */
                    if (from->id() < to->id()) {
                        if (from->time() > to->time()) {
                            swap(from, to);
                            swap(l, r);
                        }
                        int la = l->vertex()->id();
                        int lb = l->opposite()->vertex()->id();
                        int ra = r->vertex()->id();
                        int rb = r->opposite()->vertex()->id();
                        if (la == ra || la == rb || lb == ra || lb == rb) {
							this->skeleton_otherlines.push_back(Segment_2(from->point(), to->point()));
                        }
                        else {  // centerline
							this->skeleton_centerlines.push_back(Segment_2(from->point(), to->point()));
                            centerline_vertex2cnt[from] += 1;
							centerline_vertex2cnt[to] += 1;
                        }
                    }
                }
            }
			
            // find connected vertex and leaf vertex
            std::unordered_set<Vertex_handle> connected_vertex;
            std::unordered_set<Vertex_handle> leaf_vertex;
            for (auto it : centerline_vertex2cnt)
                if (it.second >= 3)
                    connected_vertex.insert(it.first);
                else if (it.second == 1)
                    leaf_vertex.insert(it.first);
			std::cout << "connected_vertex.size() = " << connected_vertex.size() << ", leaf_vertex.size() = " << leaf_vertex.size() << std::endl;

            // init partition
			auto get_adj_centerline_vertex = [&](Vertex_handle v) {
				std::vector<Vertex_handle> adj_vertex;
				Halfedge_handle h = v->halfedge();
                do {
					Vertex_handle adj = h->opposite()->vertex();
					if (centerline_vertex2cnt.find(adj) != centerline_vertex2cnt.end())
						adj_vertex.push_back(adj);
					h = h->next()->opposite();
				} while (h != v->halfedge());
				return adj_vertex;
			};
            std::vector<std::unordered_set<Vertex_handle>> parts;
            std::unordered_set<Vertex_handle> visited; 
            for (Vertex_handle it : connected_vertex) {
                std::vector<Vertex_handle> adj_vertex = get_adj_centerline_vertex(it);
                //std::cout << adj_vertex.size() << " ";
                for (auto v : adj_vertex) { 
                    if(visited.find(v) != visited.end())
						continue;
                    std::unordered_set<Vertex_handle> part;
                    part.insert(it);
                    part.insert(v);
                    visited.insert(it);
                    visited.insert(v);
                    auto adj_v = get_adj_centerline_vertex(v);
                    while (adj_v.size() == 2) {
                        auto next_v = (part.find(adj_v[0]) == part.end()) ? adj_v[0] : adj_v[1];
                        part.insert(next_v);
                        visited.insert(next_v);
						adj_v = get_adj_centerline_vertex(next_v);
                    }
                    parts.push_back(part);
                }
            }
            std::cout << "partition num : " << parts.size() << ", parition edges num : ";
            for (auto p : parts)
                std::cout << p.size() << " ";
            std::cout << std::endl;
            
            std::unordered_map < Vertex_handle, std::unordered_set<int>> v2pn;  // vertex 2 part num
            for (int part_num = 0; part_num < parts.size(); ++part_num) {
                for (auto v : parts[part_num])
                    v2pn[v].insert(part_num);
            }
			this->init_partition = std::vector<std::vector<Polygon_2>>(parts.size());
			std::vector<std::unordered_set<Face_handle>> certain_faces(parts.size());     // faces which are certain belong to which part
			std::vector<Face_handle> uncertain_faces;                       // faces which are uncertain belong to which part   
            for (Face_iterator it = this->skeleton->faces_begin(); it != this->skeleton->faces_end(); ++it) {
				Face_handle face = it;
                std::vector<Vertex_handle> centerline_vertex_in_face;
                Halfedge_handle h = face->halfedge();
                do {
                    Vertex_handle v = h->vertex();
                    if (centerline_vertex2cnt.find(v) != centerline_vertex2cnt.end())
                        centerline_vertex_in_face.push_back(v);
                    h = h->next();
                } while (h != face->halfedge());
                // all centerline is point to a center part
                int part_num = -1;
                for (auto v : centerline_vertex_in_face) {
                    if (v2pn[v].size() == 1) {
                        if(part_num == -1)
						    part_num = *v2pn[v].begin();
						else if (part_num != *v2pn[v].begin()) {
							part_num = -1;
							break;
						}
                    }
                    else {
                        part_num = -1;
                        break;
                    }
                }
                if (part_num != -1) {  // certain
					certain_faces[part_num].insert(face);
                    init_partition[part_num].push_back(get_Poly_from_Face(face));
                }
                else {  // uncertain
                    uncertain_faces.push_back(face);
					uncertain_parts.push_back(get_Poly_from_Face(face));
                }
			}

            // spilit uncertain faces
            auto areTwoIntersectionPointsInside = [](const Polygon_2& polygon, const Segment_2& segment) {
                std::vector<Point_2> intersection_points;
                Point_2 mid_point = segment.source() + (segment.target() - segment.source()) / 2;
                if (!(polygon.bounded_side(mid_point) == CGAL::ON_BOUNDED_SIDE))
                    return false;

                for (auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge) {
                    auto result = CGAL::intersection(*edge, segment);
                    if (result) {
                        if (const Point_2* p = boost::get<Point_2>(&*result)) {
                            if (std::find(intersection_points.begin(), intersection_points.end(), *p) == intersection_points.end()) {
                                intersection_points.push_back(*p);
                            }
                        }
                    }
                    if (intersection_points.size() > 2) {
                        return false;
                    }
                }
                return intersection_points.size() == 2;
                };
            for (Face_handle& face : uncertain_faces) {
                Polygon_2 poly = get_Poly_from_Face(face);
                Halfedge_handle border_edge;
                struct VAttr {  // vertex attribute
                    FT dis_to_A = 0;
                    FT dis_to_B = 0;
                    bool IsIntersect_to_A = false;
                    bool IsIntersect_to_B = false;
                    int ChooseNum() {   // 0: both not, 1: A, 2: B
                        if (IsIntersect_to_A == true && IsIntersect_to_B == true)
                            return 0;
                        else if (IsIntersect_to_B == true)
                            return 1;
                        else if (IsIntersect_to_A == true)
                            return 2;
                        else
                            return dis_to_A < dis_to_B ? 1 : 2;
                    }
                    int chooseNum;
                    Vertex_handle A, B;
                    bool operator<(const VAttr& other) const {
						return chooseNum > other.chooseNum; // 2->1->0
					}
                };
                std::unordered_map<Vertex_handle, VAttr> split_vertexs;
                Halfedge_handle h = face->halfedge();
                do {
                    Vertex_handle v = h->vertex();
                    if (connected_vertex.find(v) != connected_vertex.end())
                        split_vertexs[v] = VAttr();
                    if (!h->is_bisector()) {
                        border_edge = h;
                    }
                    h = h->next();
                } while (h != face->halfedge());
                // cal A
                Vertex_handle A = border_edge->vertex();
                Halfedge_handle border_next_halfedge = border_edge->next();
                h = border_next_halfedge;
                FT dis = 0;
                while (h != border_edge) {
                    Vertex_handle v = h->vertex();
                    dis += CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), h->prev()->vertex()->point()));
                    if (split_vertexs.find(v) != split_vertexs.end()) {
                        split_vertexs[v].dis_to_A = dis;
                        Segment_2 seg(v->point(), A->point());
                        split_vertexs[v].IsIntersect_to_A = !areTwoIntersectionPointsInside(poly, seg);
                    }
                    h = h->next();
                };
                // cal B
                Vertex_handle B = border_edge->prev()->vertex();
                Halfedge_handle border_prev_halfedge = border_edge->prev();
                h = border_prev_halfedge;
                dis = 0;
                while (h != border_edge) {
                    Vertex_handle v = h->prev()->vertex();
                    dis += CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), h->vertex()->point()));
                    if (split_vertexs.find(v) != split_vertexs.end()) {
                        split_vertexs[v].dis_to_B = dis;
                        Segment_2 seg(v->point(), B->point());
                        split_vertexs[v].IsIntersect_to_B = !areTwoIntersectionPointsInside(poly, seg);
                    }
                    h = h->prev();
                };
                // get choose num
                for (auto& it : split_vertexs) {
                    it.second.chooseNum = it.second.ChooseNum();
                    if(it.second.chooseNum == 1)
						it.second.A = A;
					else if (it.second.chooseNum == 2)
						it.second.B = B;
                }
                // sort according to choose num, B first, then A, then both not
                std::vector<std::pair<Vertex_handle, VAttr>> split_vertexs_vec(split_vertexs.begin(), split_vertexs.end());
                std::sort(split_vertexs_vec.begin(), split_vertexs_vec.end());
                for (auto& sv : split_vertexs_vec) {
                    Vertex_handle& split_v = sv.first;
                    VAttr& attr = sv.second;
                    std::cout << "chooseNum = " << attr.chooseNum << ", dis_to_A = " << attr.dis_to_A << ", dist_to_B = " << attr.dis_to_B << ", IsIntersect_to_A = " << attr.IsIntersect_to_A << ", IsIntersect_to_B = " << attr.IsIntersect_to_B << std::endl;
                    Polygon_2 poly_paint;
                    int part_num = -1;
                    if (attr.chooseNum == 1) {
                        Vertex_handle A = attr.A;
                        Halfedge_handle cur_h = A->halfedge();
                        Vertex_handle cur_v;
                        do {
                            cur_v = cur_h->vertex();
                            if (v2pn[cur_v].size() == 1)
                                part_num = *v2pn[cur_v].begin();
                            poly_paint.push_back(cur_v->point());
                            cur_h = cur_h->next();
                        } while (cur_v != split_v);
                    }
                    else if (attr.chooseNum == 2) {
                        Halfedge_handle cur_h = border_edge;
                        Vertex_handle cur_v;
                        do {
                            cur_v = cur_h->prev()->vertex();
                            cur_h = cur_h->prev();
                        } while (cur_v != split_v);
                        do {
                            cur_v = cur_h->vertex();
                            if (v2pn[cur_v].size() == 1) {
                                part_num = *v2pn[cur_v].begin();
                            }
                            poly_paint.push_back(cur_v->point());
                            cur_h = cur_h->next();
                        } while (cur_v != border_edge->prev()->vertex());
                    }
                    else {
                        continue;
                    }
                    // TODO: another part
                    if (part_num != -1)
                        init_partition[part_num].push_back(poly_paint);
                }
            }
            
            std::cout << "partition done!" << std::endl;
        }
		else
			polygon = origin_space;
    }
} // namespace Partition
#endif // PARTITION_HPP