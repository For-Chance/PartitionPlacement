#include "CenterLineContext.hpp"
#include "GeoJSON.hpp"
#include "Viewer.hpp"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"

namespace CenterLine {
    struct CenterLineTest {
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;
        using Solver = Solver<K>;
        using Context = Context<K>;

        using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
        using Polygon_2 = CGAL::Polygon_2<K>;
        using Point_2 = CGAL::Point_2<K>;
        using Vector_2 = CGAL::Vector_2<K>;

        /*
        * @brief Parameters of controlling display.
        */
        struct ControlProps {
            bool showId;	/* Show all id of points whether or not. */
            int showMode;	/* false: show polygon; true: show polygon and centerline. */
            std::string InputFile;
            std::string OutputFile;
        };

        CenterLineTest(const std::string& config_file);

        void parseConfigFile(const std::string& config_file, Context& context, ControlProps& controlProps);
        std::string parseInputFile(const std::string& input_file);
        void showResult(const Solver& clSolver);
        void write_geojson(const std::string& geojson, const std::string& outputPath);

        std::string geojson;
        Context context;
        ControlProps controlProps;
    };
}

namespace CenterLine {
    CenterLineTest::CenterLineTest(const std::string& config_file) {
        parseConfigFile(config_file, context, controlProps);
        geojson = parseInputFile(controlProps.InputFile);
        Solver clSolver(geojson, context.centerlineProps);
        showResult(clSolver);
    }

    void CenterLineTest::parseConfigFile(const std::string& config_file, Context& context, ControlProps& controlProps) {
        Json::Value p = GeoJSON::parse_propsjson(config_file);
        std::string filename = p["filename"].asString();
        controlProps.showId = p["Cprops"]["showId"].asBool();
        controlProps.showMode = p["Cprops"]["showMode"].asInt();
        controlProps.InputFile = p["Cprops"]["InputFile"].asString();
        controlProps.OutputFile = p["Cprops"]["OutputFile"].asString();
	}

    std::string CenterLineTest::parseInputFile(const std::string& input_file_path) {
        std::ifstream in(input_file_path);
        if (!in.is_open()) throw(input_file_path + "not found");
        std::istreambuf_iterator<char>  beg(in), end;
        return std::string(beg, end);
    }

    void CenterLineTest::showResult(const Solver& clSolver) {
		std::vector<Point_2> points;
		std::vector<std::pair<int, int>> segs, new_segs, sub_segs;
		std::vector<Polygon_2> PolyParts;

		const Polygon_with_holes_2& space = clSolver.origin_space; // continue
		points.clear();
		segs.clear();
		new_segs.clear();
		sub_segs.clear();
		PolyParts.clear();

		// the polygon with holes is stored in PolyParts
		PolyParts.push_back(space.outer_boundary());
		PolyParts.insert(PolyParts.end(), space.holes_begin(), space.holes_end());

		if (controlProps.OutputFile != "" && controlProps.OutputFile.substr(controlProps.OutputFile.length() - 7, 7) == "geojson") {
			write_geojson(clSolver.centerline_geojson(), controlProps.OutputFile);
		}

		if (context.centerlineProps.isSmooth == false) {
			auto result = clSolver.centerline();
			auto sub_line = clSolver.sub_centerline();

			// # get all points
			for (auto& seg : result) {
				points.push_back(seg.source());
				points.push_back(seg.target());
			}
			for (auto& seg : sub_line) {
				points.push_back(seg.source());
				points.push_back(seg.target());
			}
			std::sort(points.begin(), points.end());
			auto unique_end = std::unique(points.begin(), points.end());
			points.erase(unique_end, points.end());

			// # get segs (segments before connect), new_segs (segments after connect) and sub_segs (from sub_line)
			for (size_t i = 0; i < clSolver.seg_cnt_before_connect; ++i) {
				int a = std::lower_bound(points.begin(), points.end(), result[i].source()) - points.begin();
				int b = std::lower_bound(points.begin(), points.end(), result[i].target()) - points.begin();
				segs.emplace_back(a, b);
			}
			for (size_t i = clSolver.seg_cnt_before_connect; i < result.size(); ++i) {
				int a = std::lower_bound(points.begin(), points.end(), result[i].source()) - points.begin();
				int b = std::lower_bound(points.begin(), points.end(), result[i].target()) - points.begin();
				new_segs.emplace_back(a, b);
			}
			for (auto& seg : sub_line) {
				int a = std::lower_bound(points.begin(), points.end(), seg.source()) - points.begin();
				int b = std::lower_bound(points.begin(), points.end(), seg.target()) - points.begin();
				sub_segs.emplace_back(a, b);
			}

			// # Qt gui show
			int argc = 1;
			const char* argv[2] = { "t2_viewer","\0" };
			QApplication app(argc, const_cast<char**>(argv));
			std::string title = "output-" + controlProps.InputFile;
			CGAL::CenterLineViewer<K> mainwindow(app.activeWindow(), *PolyParts.begin(), title.c_str());
			mainwindow.drawPartitions(PolyParts, CGAL::Color(0, 0, 0), CGAL::Color(67, 177, 235));
			mainwindow.drawTree(points, segs, CGAL::Color(255, 255, 0));	// centerline
			mainwindow.drawTree(points, new_segs, CGAL::Color(255, 0, 0));	// connect line of centerline
			mainwindow.drawTree(points, sub_segs, CGAL::Color(0, 255, 0));
			mainwindow.show();
			if (controlProps.OutputFile != "") {
				std::cout << "save image" << std::endl;
				std::string outputPath = controlProps.OutputFile + title + ".svg";
				mainwindow.saveImage(outputPath);
				std::cout << "save image end: " << outputPath << std::endl;
			}
			app.exec();
		}
		else {
			auto smooth_result = clSolver.smooth_centerline();

			// # get all points
			for (auto& seg : smooth_result) {
				points.push_back(seg.source());
				points.push_back(seg.target());
			}
			std::sort(points.begin(), points.end());
			auto unique_end = std::unique(points.begin(), points.end());
			points.erase(unique_end, points.end());

			for (size_t i = 0; i < smooth_result.size(); ++i) {
				int a = std::lower_bound(points.begin(), points.end(), smooth_result[i].source()) - points.begin();
				int b = std::lower_bound(points.begin(), points.end(), smooth_result[i].target()) - points.begin();
				segs.emplace_back(a, b);
			}

			// # Qt gui show
			int argc = 1;
			const char* argv[2] = { "t2_viewer","\0" };
			QApplication app(argc, const_cast<char**>(argv));
			std::string title = "output-" + controlProps.InputFile;
			CGAL::CenterLineViewer<K> mainwindow(app.activeWindow(), *PolyParts.begin(), title.c_str());
			mainwindow.drawPartitions(PolyParts, CGAL::Color(0, 0, 0), CGAL::Color(67, 177, 235));
			mainwindow.drawTree(points, segs, CGAL::Color(255, 255, 0));	// centerline
			mainwindow.show();
			if (controlProps.OutputFile != "") {
				std::cout << "save image" << std::endl;
				std::string outputPath = controlProps.OutputFile + title + ".svg";
				mainwindow.saveImage(outputPath);
				std::cout << "save image end: " << outputPath << std::endl;
			}
			app.exec();
		}
    }

    void CenterLineTest::write_geojson(const std::string& geojson, const std::string& outputPath) {
        std::ofstream outputFile(outputPath);
        if (outputFile.is_open()) {
            outputFile << geojson;
            outputFile.close();
            std::cout << "Write completed: " << outputPath << std::endl;
        }
        else {
            throw "Cannot open the path: " + outputPath;
        }
    }
}