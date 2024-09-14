#include "Partition.hpp"
#include "Placement.hpp"
#include "PartitionPlacementContext.hpp"
#include "GeoJSON.hpp"
#include "Viewer.hpp"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Segment_2.h>

namespace ParitionPlacement
{
	struct PartitionPlacementTest
	{
		using K = CGAL::Exact_predicates_exact_constructions_kernel;
		using PartitionProps = Partition::PartitionProps<K>;
		using PlacementProps = Placement::PlacementProps<K>;

		using Polygon_2 = CGAL::Polygon_2<K>;
		using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
		using Segment_2 = CGAL::Segment_2<K>;
		using Point_2 = CGAL::Point_2<K>;

		struct ControlProps {
			std::string InputFile;
			std::string OutputFile;
			bool withSimplifyBoundary;
			bool showSkeletonFaces;
			bool showCenterLines;
			bool showPartitions;
		};

		PartitionPlacementTest(const std::string& config_file);
		void parseConfigFile(const std::string& config_file, Context<K>& context, ControlProps& controlProps);
		std::string parseInputFile(const std::string& input_file);
		void showResult(const Solver<K>& ppSolver);

		std::string geojson;
		Context<K> context;
		ControlProps controlProps;
	};
}

namespace ParitionPlacement
{
	/// <summary>
	/// Constructor: 
	///     1. parse config file; 
	///     2. parse input file; 
	///     3. solve; 
	///     4. show result
	/// </summary>
	/// <param name="config_file"></param>
	PartitionPlacementTest::PartitionPlacementTest(const std::string& config_file) {
		std::cout << "PartitionPlacementTest(const std::string& config_file)" << std::endl;
		parseConfigFile(config_file, context, controlProps);
		geojson = parseInputFile(controlProps.InputFile);
		Solver<K> ppSolver(geojson, context.ppProps.partProps, context.ppProps.placeProps);
		showResult(ppSolver);
	}

	/// <summary>
	/// Parse intput file to string
	/// </summary>
	/// <param name="input_file_path"></param>
	/// <returns></returns>
	std::string PartitionPlacementTest::parseInputFile(const std::string& input_file_path) {
		std::ifstream in(input_file_path);
		if (!in.is_open()) throw(input_file_path + "not found");
		std::istreambuf_iterator<char>  beg(in), end;
		return std::string(beg, end);
	}

	/// <summary>
	/// Parse config file to PartitionProps, PlacementProps, ControlProps
	/// </summary>
	/// <param name="config_file"></param>
	/// <param name="context"></param>
	/// <param name="controlProps"></param>
	void PartitionPlacementTest::parseConfigFile(const std::string& config_file, Context<K>& context, ControlProps& controlProps) {
		Json::Value p = GeoJSON::parse_propsjson(config_file);
		context.ppProps.partProps.withSimplifyBoundary = p["PartitionProps"]["withSimplifyBoundary"].asBool();
		context.ppProps.partProps.simplify_order = p["PartitionProps"]["simplify_order"].asString();
		controlProps.InputFile = p["ControlProps"]["InputFile"].asString();
		controlProps.OutputFile = p["ControlProps"]["OutputFile"].asString();
		controlProps.withSimplifyBoundary = p["ControlProps"]["withSimplifyBoundary"].asBool();
		controlProps.showSkeletonFaces = p["ControlProps"]["showSkeletonFaces"].asBool();
		controlProps.showCenterLines = p["ControlProps"]["showCenterLines"].asBool();
		controlProps.showPartitions = p["ControlProps"]["showPartitions"].asBool();
	}

	/// <summary>
	/// Show result with QT
	/// </summary>
	/// <param name="ppSolver"></param>
	void PartitionPlacementTest::showResult(const Solver<K>& ppSolver) {
		std::cout << "Show Result:" << std::endl;
		const Polygon_with_holes_2& space = ppSolver.origin_space;
		const std::vector<Segment_2>& skeleton_centerlines = ppSolver.get_skeleton_centerlines();
		const std::vector<Segment_2>& skeleton_otherlines = ppSolver.get_skeleton_otherlines();
		const std::vector<Polygon_2>& skeleton_faces = ppSolver.get_skeleton_faces();
		const std::vector<std::vector<Polygon_2>>& partition = ppSolver.get_partition();
		const std::vector<Polygon_2>& uncertain_parts = ppSolver.get_uncertain_parts();
		const std::vector<Segment_2>& split_segments = ppSolver.get_split_segments();
		std::vector<Polygon_2> PolyParts_outer, PolyParts_holes;
		PolyParts_outer.push_back(space.outer_boundary());
		PolyParts_holes.insert(PolyParts_holes.end(), space.holes_begin(), space.holes_end());

		// # Qt gui show
		int argc = 1;
		const char* argv[2] = { "viewer","\0" };
		QApplication app(argc, const_cast<char**>(argv));
		int screenWidth = 1920;
		int screenHeight = 1080;
		int centerX2 = screenWidth * 2 / 4;

		CGAL::Color Point_Color(0, 0, 0);
		CGAL::Color Room_Color(67, 177, 235);
		CGAL::Color Holes_Color(255, 255, 255);
		CGAL::Color Segment_Color(50, 50, 50);
		CGAL::Color Skeleton_Color(0, 255, 0);
		CGAL::Color Centerline_Color(255, 255, 0);
		CGAL::Color SplitSeg_Color(255, 0, 0);
		CGAL::Color Log_Color(0, 0, 255);

		std::string title = "output-" + controlProps.InputFile;
		CGAL::CenterLineViewer<K> mainwindow(app.activeWindow(), *PolyParts_outer.begin(), title.c_str());
		mainwindow.move(centerX2 - mainwindow.width() / 2, screenHeight / 2 - mainwindow.height() / 2);
		if (controlProps.showSkeletonFaces) {
			mainwindow.drawPartitions_withRandomColor(skeleton_faces);
		}
		if (controlProps.showCenterLines) {
			mainwindow.drawSegments(skeleton_centerlines, Centerline_Color);
			mainwindow.drawSegments(skeleton_otherlines, Skeleton_Color);
		}
		if (controlProps.showPartitions) {
			/*mainwindow.drawPartitions(PolyParts_holes, Segment_Color, Holes_Color, Point_Color);
			mainwindow.drawPartitions(PolyParts_outer, Segment_Color, Room_Color, Point_Color);*/

			for (int part_num = 0; part_num < partition.size(); part_num++) {
				const std::vector<Polygon_2>& part = partition[part_num];
				if (part.size() == 0)
					continue;
				CGAL::Color part_Color(rand() % 255, rand() % 255, rand() % 255);
				mainwindow.drawPartitions(part, Segment_Color, part_Color, Point_Color);
				// draw part_num in the center of the partition
				Point_2 center(0, 0);
				int vertex_num = 0;
				for (const Polygon_2& poly : part)
					for (const Point_2& it : poly.vertices()) {
						vertex_num++;
						center = Point_2(center.x() + it.x(), center.y() + it.y());
					}
				center = Point_2(center.x() / vertex_num, center.y() / vertex_num);
				mainwindow.drawText(center, QString::number(part_num));
			}
			mainwindow.drawPartitions(uncertain_parts, Segment_Color, Room_Color, Point_Color);
			mainwindow.drawSegments(split_segments, SplitSeg_Color);
		}
		if (controlProps.withSimplifyBoundary) {
			CGAL::Color Simplify_Segment_Color(150, 150, 150);
			const Polygon_with_holes_2& polygon = ppSolver.polygon;
			std::vector<Polygon_2> PolyParts;
			PolyParts.push_back(polygon.outer_boundary());
			PolyParts.insert(PolyParts.end(), polygon.holes_begin(), polygon.holes_end());
			std::vector<Segment_2> segs;
			for (auto& poly : PolyParts)
				for (auto it : poly.edges())
					segs.push_back(it);
			mainwindow.drawSegments(segs, Simplify_Segment_Color);
		}
		mainwindow.show();
		app.exec();
	}
}