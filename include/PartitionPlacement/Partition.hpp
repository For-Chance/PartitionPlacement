#ifndef PARTITION_HPP
#define PARTITION_HPP
#include "stdafx.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> 
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <list>
#include <chrono>

#include "SimplifyBoundary.hpp"
#include "convertKernel.hpp"

namespace Partition
{
	const double eps = 1. / (1 << 20), eps2 = eps * eps,
		eps_limit = 0.25 / (1ll << 62);

	template <typename K>
	struct PartitionProps {
		bool withSimplifyBoundary;
		std::string simplify_order;
		PartitionProps() {
			withSimplifyBoundary = true;
			simplify_order = "0";
		}
	};


	template <typename K>
	class Solver
	{
		using InnerK = CGAL::Exact_predicates_inexact_constructions_kernel;
		using FT = typename InnerK::FT;
		using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<InnerK>;
		using Polygon_2 = CGAL::Polygon_2<InnerK>;
		using Point_2 = CGAL::Point_2<InnerK>;
		using Segment_2 = CGAL::Segment_2<InnerK>;
		using Polygon_list = std::list<Polygon_2>;
		using Vector_2 = CGAL::Vector_2<InnerK>;
		using Line_2 = CGAL::Line_2<InnerK>;
		using Polygon_set_2 = CGAL::Polygon_set_2<InnerK>;
		using Ss = CGAL::Straight_skeleton_2<InnerK>;
		using Halfedge_iterator = typename Ss::Halfedge_iterator;
		using Halfedge_handle = typename Ss::Halfedge_handle;
		using Vertex_handle = typename Ss::Vertex_handle;
		using Vertex_iterator = typename Ss::Vertex_iterator;
		using Face_handle = typename Ss::Face_handle;
		using Face_iterator = typename Ss::Face_iterator;
		using HDS = typename Ss::Base;
		using Decorator = CGAL::HalfedgeDS_decorator<HDS>;
		using SsBuilderTraits = CGAL::Straight_skeleton_builder_traits_2<InnerK>;
		using SsBuilder = CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss>;

		using PartitionProps = PartitionProps<K>;
	private:
		struct Point2Hash {
			std::size_t operator()(const Point_2& p) const {
				std::size_t seed = 0;
				boost::hash<double> hasher;
				auto add_hash = [&seed, &hasher](double value) {
					seed ^= hasher(value) + 0x9e3129b9 + (seed << 6) + (seed >> 2);
					};

				add_hash(CGAL::to_double(p.x()));
				add_hash(CGAL::to_double(p.y()));

				return seed;
			}
		};
		struct Point2Equal {
			bool operator()(const Point_2& p1, const Point_2& p2) const {
				return p1.x() == p2.x() && p1.y() == p2.y();
			}
		};
		static bool equal(const FT& a, const FT& b, FT epsilon = eps)
		{
			return a + epsilon > b && a < b + epsilon;
		}
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
			std::cout << "skeleton built! Time elapsed:" << diff.count() << " s, num of skeleton face:" << skeleton->size_of_faces() << std::endl;
			return skeleton;
		};
		/// <summary>
		/// get  Polygon_2 from a Face_handle
		/// </summary>
		/// <param name="face"></param>
		/// <returns></returns>
		Polygon_2 get_Poly_from_Face(const Face_handle& face) {
			std::vector<Point_2> ps;
			ps.reserve(20);
			Halfedge_handle h = face->halfedge();
			Halfedge_handle first = h;
			do {
				Point_2 point = h->vertex()->point();
				ps.push_back(point);
				h = h->next();
			} while (h != first);
			return Polygon_2(ps.begin(), ps.end());
		}

		KernelConverter::KernelConverter<K, InnerK, KernelConverter::NumberConverter<typename K::FT, FT>>K2InnerK;
		KernelConverter::KernelConverter<InnerK, K, KernelConverter::NumberConverter<FT, typename K::FT>>InnerK2K;
		CGAL::Polygon_with_holes_2<K> origin_space;
		Polygon_with_holes_2 polygon;
		boost::shared_ptr<Ss> skeleton;
		std::vector<Segment_2> skeleton_segments;
		std::vector<Segment_2> skeleton_centerlines;
		std::vector<Segment_2> skeleton_otherlines;
		std::vector<Polygon_2> skeleton_faces;
		std::vector<std::vector<Polygon_2>> partition;
		std::vector<Polygon_2> merge_partition;
		std::vector<Polygon_2> uncertain_parts;
		std::vector<Segment_2> split_segments;
		std::vector<Point_2> log_points;
	public:
		Solver() {}
		Solver(const CGAL::Polygon_with_holes_2<K>& space, const PartitionProps& props = PartitionProps());

		CGAL::Polygon_with_holes_2<K> get_origin_space() const { return this->origin_space; }
		CGAL::Polygon_with_holes_2<K> get_polygon() const { return this->InnerK2K.convert(this->polygon); }
		std::vector<CGAL::Segment_2<K>> get_skeleton_segments() const { this->InnerK2K.convert(this->skeleton_segments); }
		std::vector<CGAL::Segment_2<K>> get_skeleton_centerlines() const { return this->InnerK2K.convert(this->skeleton_centerlines); }
		std::vector<CGAL::Segment_2<K>> get_skeleton_otherlines() const { return this->InnerK2K.convert(this->skeleton_otherlines); }
		std::vector<CGAL::Polygon_2<K>> get_skeleton_faces() const { return this->InnerK2K.convert(this->skeleton_faces); }
		std::vector<std::vector<CGAL::Polygon_2<K>>> get_partition() const {
			std::vector<std::vector<CGAL::Polygon_2<K>>> res;
			for (const auto& part : this->partition)
				res.push_back(this->InnerK2K.convert(part));
			return res;
		}
		std::vector<CGAL::Polygon_2<K>> get_uncertain_parts() const { return this->InnerK2K.convert(this->uncertain_parts); }
		std::vector<CGAL::Segment_2<K>> get_split_segments() const { return this->InnerK2K.convert(this->split_segments); }
		std::vector<CGAL::Point_2<K>> get_log_points() const { return this->InnerK2K.convert(this->log_points); }
	};
}

namespace Partition
{
	template <typename K>
	Solver<K>::Solver(const CGAL::Polygon_with_holes_2<K>& space, const PartitionProps& props)
	{
		this->origin_space = space;
		if (props.withSimplifyBoundary) {
			// 1. Simplify Boundary
			SimplifyBoundary::SimplifyProps<InnerK> simpProps = SimplifyBoundary::SimplifyProps<InnerK>();
			SimplifyBoundary::ExpandProps<InnerK> expandProps = SimplifyBoundary::ExpandProps<InnerK>();
			expandProps.simplify_order = props.simplify_order;
			SimplifyBoundary::Solver<InnerK> sbSolver(this->K2InnerK.convert(this->origin_space), simpProps, expandProps);
			this->polygon = sbSolver.simplify_space;

			// 2. skeleton 
			this->skeleton = build_skeleton(this->polygon);
			for (Halfedge_iterator it = this->skeleton->halfedges_begin(); it != this->skeleton->halfedges_end(); ++it) {
				Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
				skeleton_segments.push_back(Segment_2(from->point(), to->point()));
			}
			Decorator decorator(*this->skeleton);

			std::cout << "partition..." << std::endl;
			// 3. find centerline
			// 3.1 save to skeleton_faces
			for (Face_iterator it = this->skeleton->faces_begin(); it != this->skeleton->faces_end(); ++it) {
				Face_handle face = it;
				skeleton_faces.push_back(get_Poly_from_Face(face));
			}
			// 3.2 find center line
			// first we try to make skeleton of no connect boundary is center line
			std::unordered_map<Vertex_handle, int> centerline_vertex2cnt;
			std::unordered_set<Vertex_handle> stop_vertices;
			std::unordered_set<Halfedge_handle> centerline_hes;
			std::vector<std::pair<Vertex_handle, bool>> he_pool;   // from_vertex 2 is_ans ; use a vector to control the loop sequence so that there is always same result
			double THRE_THETA = 150;
			FT MIN_SQUARED_COS_THETA = std::cos(THRE_THETA * M_PI / 180);
			MIN_SQUARED_COS_THETA *= MIN_SQUARED_COS_THETA;
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
						/*int la = l->vertex()->id();
						int lb = l->opposite()->vertex()->id();
						int ra = r->vertex()->id();
						int rb = r->opposite()->vertex()->id();
						if (la == ra || la == rb || lb == ra || lb == rb) {
							this->skeleton_otherlines.push_back(Segment_2(from->point(), to->point()));
						}*/
						Point_2 v_la = l->vertex()->point();
						Point_2 v_lb = l->opposite()->vertex()->point();
						Vector_2 v_l = v_la - v_lb;
						Point_2 v_ra = r->vertex()->point();
						Point_2 v_rb = r->opposite()->vertex()->point();
						Vector_2 v_r = v_ra - v_rb;
						FT inner = v_l.x() * v_r.x() + v_l.y() * v_r.y();
						FT outer = v_l.x() * v_r.y() - v_l.y() * v_r.x();
						FT absolute_scale = v_l.squared_length() * v_r.squared_length();
						FT squared_cos_theta = (inner * inner) / absolute_scale;
						bool is_ans = !(equal(v_la.x(), v_rb.x()) && equal(v_la.y(), v_rb.y()) && squared_cos_theta >= MIN_SQUARED_COS_THETA && inner >= 0);
						if (is_ans)
						{
							centerline_vertex2cnt[from] += 1;
							centerline_vertex2cnt[to] += 1;
							bool is_stop = inner < 0 && squared_cos_theta >= MIN_SQUARED_COS_THETA;
							if (is_stop) {
								stop_vertices.insert(from);
								stop_vertices.insert(to);
							}
						}
						he_pool.push_back(std::make_pair(from, is_ans));
					}
				}
			}

			// 3.3 find connected vertex and leaf vertex
			auto get_adj_centerline_vertex = [&](Vertex_handle v) {
				std::vector<Vertex_handle> adj_vertex;
				Halfedge_handle h = v->halfedge();
				do {
					Vertex_handle adj = h->opposite()->vertex();
					if (h->is_bisector() && centerline_vertex2cnt.find(adj) != centerline_vertex2cnt.end()) {
						adj_vertex.push_back(adj);
					}
					h = h->next()->opposite();
				} while (h != v->halfedge());
				return adj_vertex;
				};
			auto get_adj_centerline_halfedge = [&](Vertex_handle v) {
				std::vector<Halfedge_handle> adj_hes;
				Halfedge_handle h = v->halfedge();
				do {
					Halfedge_handle adj_h = h->opposite();
					if (h->is_bisector() && centerline_vertex2cnt.find(adj_h->vertex()) != centerline_vertex2cnt.end())
						adj_hes.push_back(adj_h);
					h = h->next()->opposite();
				} while (h != v->halfedge());
				return adj_hes;
				};
			auto get_ans_num = [&](Vertex_handle v) {
				return get_adj_centerline_vertex(v).size();
				};
			auto get_next_ans_vertex = [&](Vertex_handle v, bool& delete_v) {
				std::vector<Vertex_handle>& adj_vertices = get_adj_centerline_vertex(v);
				if (adj_vertices.size() != 1) {
					delete_v = false;
					return v;
				}
				return adj_vertices[0];
				};
			std::unordered_set<Vertex_handle> leaf_vertex;
			for (std::pair<Vertex_handle, bool>& it : he_pool)
				if (it.second == true)   // is centerline candidate
					if (centerline_vertex2cnt[it.first] == 1)
						leaf_vertex.insert(it.first);
			for (Vertex_handle lv : leaf_vertex) {
				Vertex_handle& cur_v = lv;
				log_points.push_back(cur_v->point());
				bool delete_v = true;
				do {
					Vertex_handle next_v = get_next_ans_vertex(cur_v, delete_v);
					if (!delete_v)
						break;
					centerline_vertex2cnt.erase(cur_v);
					cur_v = next_v;
				} while (stop_vertices.find(cur_v) == stop_vertices.end());
			}
			std::unordered_set<Vertex_handle> connected_vertex;
			for (auto& it : centerline_vertex2cnt)
				it.second = get_adj_centerline_vertex(it.first).size();
			leaf_vertex.clear();
			for (auto& it : centerline_vertex2cnt)
				if (it.second == 1)
					leaf_vertex.insert(it.first);
				else if (it.second == 3)
					connected_vertex.insert(it.first);
			for (Halfedge_iterator it = this->skeleton->halfedges_begin(); it != this->skeleton->halfedges_end(); ++it) {
				if (it->is_bisector()) {    // the halfedge is not border
					Vertex_handle from = it->opposite()->vertex(), to = it->vertex();
					if (from->id() < to->id()) {
						if (from->time() > to->time())
							swap(from, to);
						if (centerline_vertex2cnt.find(from) != centerline_vertex2cnt.end() && centerline_vertex2cnt.find(to) != centerline_vertex2cnt.end()) {   // centerline
							this->skeleton_centerlines.push_back(Segment_2(from->point(), to->point()));
							centerline_hes.insert(it);
							centerline_hes.insert(it->opposite());
						}
						else
							this->skeleton_otherlines.push_back(Segment_2(from->point(), to->point()));
					}
				}
			}

			// 4. partition
			// 4.1 init partition according to connected vertex
			std::vector<std::unordered_set<Face_handle>> Face_partition;
			std::vector<std::unordered_set<Vertex_handle>> parts;
			std::unordered_set<Halfedge_handle> special_hes;
			std::unordered_set<Halfedge_handle> visited;
			for (Vertex_handle it : connected_vertex) {
				std::vector<Halfedge_handle>& adj_hes = get_adj_centerline_halfedge(it);
				for (auto he : adj_hes) {
					if (visited.find(he) != visited.end())
						continue;
					Vertex_handle v = he->vertex();
					std::unordered_set<Vertex_handle> part;
					int part_num = parts.size();
					visited.insert(he);
					visited.insert(he->opposite());
					part.insert(it);
					part.insert(v);
					auto adj_next_hes = get_adj_centerline_halfedge(v);
					while (adj_next_hes.size() == 2) {
						Halfedge_handle next_hes = (part.find(adj_next_hes[0]->vertex()) == part.end()) ? adj_next_hes[0] : adj_next_hes[1];
						if (visited.find(next_hes) != visited.end())
							break;
						visited.insert(next_hes);
						visited.insert(next_hes->opposite());
						Vertex_handle next_v = next_hes->vertex();
						part.insert(next_v);
						adj_next_hes = get_adj_centerline_halfedge(next_v);
					}
					if (part.size() == 2 && connected_vertex.find(v) != connected_vertex.end()) {
						special_hes.insert(he);
						special_hes.insert(he->opposite());
					}
					parts.push_back(part);
				}
			}
			std::cout << "partition num : " << parts.size() << ", parition edges num : ";
			int total_partition_edges_num = 0;
			for (auto p : parts) {
				std::cout << p.size() << " ";
				total_partition_edges_num += p.size();
			}
			std::cout << ", sum up: " << total_partition_edges_num << std::endl;
			std::unordered_map < Vertex_handle, std::unordered_set<int>> v2pn;  // vertex 2 part num
			for (int part_num = 0; part_num < parts.size(); ++part_num)
				for (auto v : parts[part_num])
					v2pn[v].insert(part_num);
			Face_partition = std::vector<std::unordered_set<Face_handle>>(parts.size());

			// 4.2 find certain faces and uncertain faces
			auto get_adj_faces = [&](const Face_handle& face) {
				std::unordered_set<Face_handle> adj_faces;
				Halfedge_handle h = face->halfedge();
				do {
					if (h->is_bisector())
						adj_faces.insert(h->opposite()->face());
					h = h->next();
				} while (h != face->halfedge());
				return adj_faces;
				};
			std::vector<std::unordered_set<Face_handle>> certain_faces(parts.size());     // faces which are certain belong to which part
			std::unordered_set<Face_handle> uncertain_faces;                       // faces which are uncertain belong to which part   
			std::unordered_map<Face_handle, int> face2pn;                   // face 2 part num
			std::unordered_map<Halfedge_handle, int> he2pn;
			for (int i = 0; i < parts.size(); i++) {
				const int& part_num = i;
				const std::unordered_set<Vertex_handle>& part = parts[i];
				for (auto& v : part) {
					Halfedge_handle h = v->halfedge();
					do {
						Vertex_handle from = h->opposite()->vertex(), to = h->vertex();
						if (part.find(from) != part.end() && part.find(to) != part.end() && he2pn.find(h) == he2pn.end()) {
							he2pn[h] = part_num;
							he2pn[h->opposite()] = part_num;
						}
						h = h->next()->opposite();
					} while (h != v->halfedge());
				}
			}
			for (auto it : he2pn) {
				Halfedge_handle he = it.first;
				int part_num = it.second;

				std::queue<Face_handle> q;
				std::unordered_set<Face_handle> q_set;  // make sure the face in the queue is unique, otherwise it will be endless loop
				Face_handle f = he->face();
				if (face2pn.find(f) != face2pn.end())
					continue;

				// the face is uncertain if there is a connected_vertex in there
				bool is_uncertain = false;
				Halfedge_handle h = f->halfedge();
				do {
					Vertex_handle v = h->vertex();
					if (connected_vertex.find(v) != connected_vertex.end()) {
						is_uncertain = true;
						break;
					}
					h = h->next();
				} while (h != f->halfedge());
				if (is_uncertain)
					uncertain_faces.insert(f);
				else {
					q.push(f);
					q_set.insert(f);
				}

				// the face is certain if there is a leaf vertex
				if (leaf_vertex.find(he->vertex()) != leaf_vertex.end()) {
					Face_handle f = he->next()->opposite()->face();
					q.push(f);
					q_set.insert(f);
				}

				while (!q.empty()) {
					Face_handle cur_face = q.front();
					q.pop();
					q_set.erase(cur_face);
					certain_faces[part_num].insert(cur_face);
					face2pn[cur_face] = part_num;
					Face_partition[part_num].insert(cur_face);
					std::unordered_set<Face_handle>& adj_faces = get_adj_faces(cur_face);
					//std::cout << "adj_faces.size() = " << adj_faces.size() << std::endl;
					for (auto& adj_face : adj_faces) {
						bool is_same_part = true;
						Halfedge_handle h = adj_face->halfedge();
						do {
							Vertex_handle v = h->vertex();
							if (connected_vertex.find(v) != connected_vertex.end()) {
								is_same_part = false;
								break;
							}
							h = h->next();
						} while (h != adj_face->halfedge());
						if (is_same_part && face2pn.find(adj_face) == face2pn.end() && q_set.find(adj_face) == q_set.end()) {
							q.push(adj_face);
							q_set.insert(adj_face);
						}
					}
				}
			}

			// 4.3 get halfedge 2 part num
			for (int part_num = 0; part_num < Face_partition.size(); part_num++)
				for (const Face_handle& face : Face_partition[part_num]) {
					Halfedge_handle cur_h = face->halfedge();
					do {
						he2pn[cur_h] = part_num;
						cur_h = cur_h->next();
					} while (cur_h != face->halfedge());
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
			std::unordered_set<Face_handle> temp_convert_faces;
			for (Face_handle face : uncertain_faces) {
				Polygon_2 poly = get_Poly_from_Face(face);
				Halfedge_handle border_edge;
				struct VAttr {  // vertex attribute
					FT dis_to_A = 0;
					int node_to_A = 0;
					FT dis_to_B = 0;
					int node_to_B = 0;
					bool IsIntersect_to_A = false;
					bool IsIntersect_to_B = false;
					int ChooseNum() {   // 0: both not, 1: A, 2: B
						if (IsIntersect_to_A == true && IsIntersect_to_B == true)
							return dis_to_A < dis_to_B ? -1 : -2;
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
						if (chooseNum != other.chooseNum)
							return chooseNum > other.chooseNum; // 2->1->-1->-2
						// ChooseNum == other.ChooseNum, make shortest distance first
						if (chooseNum == 1 || chooseNum == -1)
							return dis_to_A != other.dis_to_A ? dis_to_A < other.dis_to_A : dis_to_B < other.dis_to_B;
						else if (chooseNum == 2 || chooseNum == -2)
							return dis_to_B != other.dis_to_B ? dis_to_B < other.dis_to_B : dis_to_A < other.dis_to_A;
						else
							throw std::runtime_error("ERROR! Choose nothing!");
					}
				};
				std::unordered_map<Halfedge_handle, VAttr> split_hes;	// he->vertex() is split vertex
				std::unordered_set<int> split_hes_pns;
				Halfedge_handle h = face->halfedge();
				do {
					if (he2pn.find(h) != he2pn.end())
						split_hes_pns.insert(he2pn[h]);
					Vertex_handle v = h->vertex();
					if (connected_vertex.find(v) != connected_vertex.end())
						split_hes[h] = VAttr();
					if (!h->is_bisector()) {
						border_edge = h;
					}
					h = h->next();
				} while (h != face->halfedge());
				std::cout << "split_hes.size() = " << split_hes.size() << std::endl;
				// cal A
				if (border_edge == nullptr)
					continue;
				Vertex_handle A = border_edge->vertex();
				Halfedge_handle border_next_halfedge = border_edge->next();
				h = border_next_halfedge;
				FT dis = 0;
				int node_num = 0;
				while (h != border_edge) {
					Vertex_handle v = h->vertex();
					dis += CGAL::approximate_sqrt(CGAL::squared_distance(v->point(), h->prev()->vertex()->point()));
					node_num++;
					if (split_hes.find(h) != split_hes.end()) {
						split_hes[h].dis_to_A = dis;
						split_hes[h].node_to_A = node_num;
						Segment_2 seg(v->point(), A->point());
						split_hes[h].IsIntersect_to_A = !areTwoIntersectionPointsInside(poly, seg);
					}
					h = h->next();
				};
				// cal B
				Vertex_handle B = border_edge->prev()->vertex();
				Halfedge_handle border_prev_halfedge = border_edge->prev();
				h = border_prev_halfedge;
				dis = 0;
				node_num = 0;
				while (h != border_edge) {
					Vertex_handle v = h->vertex();
					if (split_hes.find(h) != split_hes.end()) {
						split_hes[h].dis_to_B = dis;
						split_hes[h].node_to_B = node_num;
						Segment_2 seg(v->point(), B->point());
						split_hes[h].IsIntersect_to_B = !areTwoIntersectionPointsInside(poly, seg);
					}
					dis += CGAL::approximate_sqrt(CGAL::squared_distance(h->prev()->vertex()->point(), v->point()));
					node_num++;
					h = h->prev();
				};
				// get choose num
				for (auto& it : split_hes) {
					it.second.chooseNum = it.second.ChooseNum();
					if (it.second.chooseNum == 1 || it.second.chooseNum == -1)
						it.second.A = A;
					else if (it.second.chooseNum == 2 || it.second.chooseNum == -2)
						it.second.B = B;
				}
				// sort according to choose num, B first, then A, then both not
				std::vector<std::pair<Halfedge_handle, VAttr>> split_hes_vec(split_hes.begin(), split_hes.end());
				std::sort(split_hes_vec.begin(), split_hes_vec.end(), [](const std::pair<Halfedge_handle, VAttr>& va1, const std::pair<Halfedge_handle, VAttr>& va2) {
					return va1.second < va2.second;
					});
				int another_part_num = -1;
				Face_handle spare_face = face;
				for (auto& sv : split_hes_vec) {
					Halfedge_handle split_he = sv.first;
					int shp_size = split_hes_pns.size();
					Vertex_handle split_v = split_he->vertex();
					//std::cout << "split_v: " << split_v->point() << std::endl;
					VAttr& attr = sv.second;
					//std::cout << "chooseNum = " << attr.chooseNum << ", dis_to_A = " << attr.dis_to_A << ", node_to_A = " << attr.node_to_A << ", dist_to_B = " << attr.dis_to_B << ", node_to_B = " << attr.node_to_B << ", IsIntersect_to_A = " << attr.IsIntersect_to_A << ", IsIntersect_to_B = " << attr.IsIntersect_to_B << ", split_v = " << split_v->point() << ", border_edge point = " << border_edge->vertex()->point() << std::endl;
					std::vector<Face_handle> face_paints;
					std::vector<std::pair<Halfedge_handle, Halfedge_handle>> split_segs;
					int part_num = -1;

					auto check_loop_next = [&](Halfedge_handle& start_he, const Line_2& supp_line) {
						// skip the point on supporting_line, find the first point not on the line
						Halfedge_handle cur_prev_h = start_he->prev();
						while (supp_line.has_on_boundary(cur_prev_h->vertex()->point())) {
							cur_prev_h = cur_prev_h->prev();
						};
						bool RIGHT_DIRECTION = supp_line.has_on_positive_side(cur_prev_h->vertex()->point());

						Halfedge_handle last_right_he = start_he;
						Halfedge_handle cur_h = start_he->next();
						Vertex_handle cur_v = cur_h->vertex();    // A->next
						while (cur_h != split_he) {
							if (supp_line.has_on_positive_side(cur_v->point()) == RIGHT_DIRECTION || supp_line.has_on_boundary(cur_v->point())) {
								if (last_right_he->next() != cur_h)
									split_segs.push_back(std::make_pair(cur_h, last_right_he));
								last_right_he = cur_h;
							}
							cur_h = cur_h->next();
							cur_v = cur_h->vertex();
						};
						if (last_right_he->next() != cur_h)
							split_segs.push_back(std::make_pair(cur_h, last_right_he));
						return last_right_he != start_he;	// true is intersection, false is no intersection
						};
					auto check_loop_prev = [&](Halfedge_handle& start_he, const Line_2& supp_line) {
						Halfedge_handle cur_next_h = start_he->next();
						while (supp_line.has_on_boundary(cur_next_h->vertex()->point())) {
							cur_next_h = cur_next_h->next();
						};
						bool RIGHT_DIRECTION = supp_line.has_on_positive_side(cur_next_h->vertex()->point());

						Halfedge_handle last_right_he = start_he;
						Halfedge_handle cur_h = last_right_he->prev();
						Vertex_handle cur_v = cur_h->vertex();    // B->prev
						while (cur_h != split_he) {
							if (supp_line.has_on_positive_side(cur_v->point()) == RIGHT_DIRECTION || supp_line.has_on_boundary(cur_v->point())) {
								if (last_right_he->prev() != cur_h)
									split_segs.push_back(std::make_pair(cur_h, last_right_he));
								last_right_he = cur_h;
							}
							cur_h = cur_h->prev();
							cur_v = cur_h->vertex();
						};
						if (last_right_he != start_he && last_right_he->prev() != cur_h)
							split_segs.push_back(std::make_pair(cur_h, last_right_he));
						return last_right_he != start_he;
						};

					if (attr.chooseNum == 1) {  // choose A
						part_num = he2pn[split_he];
						split_segs.push_back(std::make_pair(split_he, border_edge));
					}
					else if (attr.chooseNum == 2) {  // choose B
						part_num = he2pn[split_he->next()];
						split_segs.push_back(std::make_pair(border_edge->prev(), split_he));
					}
					else if (attr.chooseNum == -1) {    // choose A but there is intersection
						part_num = he2pn[split_he];
						Vertex_handle A = attr.A;
						Segment_2 A2sv = Segment_2(A->point(), split_v->point());
						Line_2 supp_line = A2sv.supporting_line();

						Segment_2 border_edge_seg(A->point(), B->point());
						//std::cout << "sv:" << split_v->point() << ", A:" << A->point() << ", B:" << B->point() << ", border_edge length^2:" << border_edge_seg.squared_length() << std::endl;
						Halfedge_handle start_he = border_edge;
						bool is_loop_next_intersection = check_loop_next(start_he, supp_line);
						if (!is_loop_next_intersection)
							bool is_loop_prev_intersection = check_loop_prev(start_he, supp_line);
					}
					else if (attr.chooseNum == -2) {    // choose B but there is intersection
						part_num = he2pn[split_he->next()];
						Vertex_handle B = attr.B;
						Segment_2 sv2B = Segment_2(split_v->point(), B->point());
						Line_2 supp_line = sv2B.supporting_line();

						Halfedge_handle start_he = border_edge->prev();
						bool is_loop_prev_intersection = check_loop_prev(start_he, supp_line);
						if (!is_loop_prev_intersection)
							bool is_loop_next_intersection = check_loop_next(start_he, supp_line);
					}
					else
						throw std::runtime_error("ERROR! Choose nothing!");
					for (auto it : split_segs) {
						Halfedge_handle split_face_he = decorator.split_face(it.first, it.second);
						Face_handle fp = split_face_he->face();
						spare_face = split_face_he->opposite()->face();
						temp_convert_faces.insert(face);
						split_hes_pns.erase(part_num);
						Face_partition[part_num].insert(fp);
						Halfedge_handle cur_h = fp->halfedge();
						do {
							he2pn[cur_h] = part_num;
							cur_h = cur_h->next();
						} while (cur_h != fp->halfedge());
					}
				}
				Halfedge_handle cur_h = spare_face->halfedge();
				int part_num = -1;
				if (split_hes_pns.size() == 1) {
					part_num = *split_hes_pns.begin();
					Face_partition[part_num].insert(spare_face);
					Halfedge_handle cur_h = spare_face->halfedge();
					do {
						he2pn[cur_h] = part_num;
						cur_h = cur_h->next();
					} while (cur_h != spare_face->halfedge());
				}
				else
					std::cout << "part_num == -1, chooseNum == spare_face" << std::endl;
			}
			for (Face_handle face : temp_convert_faces)
				uncertain_faces.erase(face);
			temp_convert_faces.clear();
			//std::cout << "special_hes.size() = " << special_hes.size() << std::endl;
			for (Halfedge_handle he : special_hes) {
				int part_num = he2pn[he];
				Face_handle face = he->face();
				Face_partition[part_num].insert(face);
				Halfedge_handle cur_h = face->halfedge();
				do {
					he2pn[cur_h] = part_num;
					cur_h = cur_h->next();
				} while (cur_h != face->halfedge());
				if (uncertain_faces.find(face) != uncertain_faces.end())
					uncertain_faces.erase(face);
			}
			this->partition = std::vector<std::vector<Polygon_2>>(Face_partition.size());
			for (int part_num = 0; part_num < Face_partition.size(); part_num++)
				for (const Face_handle& face : Face_partition[part_num])
					if (face != nullptr)
						this->partition[part_num].push_back(get_Poly_from_Face(face));
			std::cout << "Partition done! uncerteain_faces.size() = " << uncertain_faces.size() << std::endl;

			// 1. try merge all faces to one face so that we can get the big Polygon_2 of a part
			// 2. make sure whether is standard part for a part
			struct Part {
				bool is_standard = true;
				Polygon_2 border_polygon;
				std::unordered_set<Halfedge_handle> border_hes;
				std::unordered_set<Halfedge_handle> centerline_hes;
				std::unordered_set<int> adjacent_nonstandard_part; // neibor nonstandard part num, include itself
				void check_standard(const std::unordered_set<Halfedge_handle>& outer_centerline_hes, const std::unordered_map<Halfedge_handle, int>& he2pn) {
					for (const Halfedge_handle& he : border_hes)
						if (outer_centerline_hes.find(he) != outer_centerline_hes.end()) {
							is_standard = false;
							centerline_hes.insert(he);
							adjacent_nonstandard_part.insert(he2pn.at(he));
							adjacent_nonstandard_part.insert(he2pn.at(he->opposite()));
						}
				}
			};
			std::unordered_map<int, Part> pn2nonstandard_part;
			std::cout << "nonstandard parts(face size):" << std::endl;
			for (int part_num = 0; part_num < Face_partition.size(); part_num++) {
				Part part;
				std::unordered_set<Face_handle>& faces = Face_partition[part_num];
				if (faces.size() == 0) {
					std::cout << "empty part!, part_num == " << part_num << std::endl;
					continue;
				}
				Halfedge_handle cur_h;
				auto is_border_he = [&](const Halfedge_handle& cur_h) {
					return !(cur_h->is_bisector() && he2pn[cur_h] == he2pn[cur_h->opposite()]);
					};
				auto merge_faces = [&](const std::unordered_set<Face_handle>& faces, Polygon_2& result_polygon, std::unordered_set<Halfedge_handle>& result_hes) {
					result_hes.clear();
					for (const Face_handle& face : faces) {
						cur_h = face->halfedge();
						do {
							if (is_border_he(cur_h))
								break;
							cur_h = cur_h->next();
						} while (cur_h != face->halfedge());
						if (is_border_he(cur_h))
							break;
					}
					std::vector<Point_2> points;
					Halfedge_handle start_h = cur_h;
					do {
						result_hes.insert(cur_h);
						points.push_back(cur_h->vertex()->point());
						cur_h = cur_h->next();
						while (!is_border_he(cur_h)) {
							cur_h = cur_h->opposite()->next();
						};
					} while (cur_h != start_h);
					result_polygon = Polygon_2(points.begin(), points.end());
					};
				merge_faces(faces, part.border_polygon, part.border_hes);
				part.check_standard(centerline_hes, he2pn);
				if (!part.is_standard) {
					pn2nonstandard_part[part_num] = part;
					std::cout << part_num << "(face_size = " << faces.size() << ", centerline_size = " << part.centerline_hes.size() << ")" << std::endl;
				}
				this->partition[part_num] = { part.border_polygon };   // override the partition
			}
			// 3. merge all non standart parts to standar part
			// find the set of neibor nonstandard parts 
			// TODO bing cha ji
			std::vector<std::unordered_set<int>> nonstandard_parts_set; // neibor nonstandard set of part num
			std::unordered_map<int, int> part_num2nonstandard_set_num; // part num to nonstandard set num
			std::cout << "non standard parts:" << std::endl;
			for (auto it : pn2nonstandard_part) {
				const int& part_num = it.first;
				const Part& part = it.second;

				bool is_new_set = true;
				int set_num = -1;
				for (int pn : part.adjacent_nonstandard_part) {
					if (part_num2nonstandard_set_num.find(pn) != part_num2nonstandard_set_num.end()) {
						is_new_set = false;
						set_num = part_num2nonstandard_set_num[pn];
						break;
					}
				}
				if (is_new_set) {
					set_num = nonstandard_parts_set.size();
					nonstandard_parts_set.push_back(part.adjacent_nonstandard_part);
					for (int pn : part.adjacent_nonstandard_part)
						part_num2nonstandard_set_num[pn] = set_num;
				}
				else {
					for (int pn : part.adjacent_nonstandard_part) {
						nonstandard_parts_set[set_num].insert(pn);
						part_num2nonstandard_set_num[pn] = set_num;
					}
				}
			}
			for (auto parts : nonstandard_parts_set) {
				std::cout << "(";
				for (int part_num : parts)
					std::cout << part_num << " ";
				std::cout << ")" << std::endl;
			}

			// merge
			for (std::unordered_set<int> parts : nonstandard_parts_set) {
				int the_first_part_num = *parts.begin();
				std::vector<Polygon_2> polygons;
				for (int part_num : parts) {
					polygons.insert(polygons.end(), this->partition[part_num].begin(), this->partition[part_num].end());
					this->partition[part_num].clear();
				}
				this->partition[the_first_part_num] = polygons;
			}
		}
		else
			polygon = this->K2InnerK.convert(origin_space);
	}
} // namespace Partition
#endif // PARTITION_HPP