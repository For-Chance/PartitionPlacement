#include "SimplifyBoundaryContext.hpp"
#include "GeoJSON.hpp"
#include "Viewer.hpp"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"

namespace SimplifyBoundary {
    struct SimplifyBoundaryTest {
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

        SimplifyBoundaryTest(const std::string& config_file);

        void parseConfigFile(const std::string& config_file, Context& context, ControlProps& controlProps);
        std::string parseInputFile(const std::string& input_file);
        void showResult(const Solver& sbSolver);
        void write_geojson(const std::string& geojson, const std::string& outputPath);
        /**
        * @brief Get all points of Polygon_with_holes_2.
        *
        * @param space	A Polygon_with_holes_2 subject.
        *
        * @return A vecotr of points with its id.
        **/
        static std::vector<std::pair<Point_2, int>> get_point_id(const Polygon_with_holes_2& space) {
            std::vector<std::pair<Point_2, int>> pointIdVec;
            std::vector<Polygon_2> PolyParts;
            PolyParts.push_back(space.outer_boundary());
            PolyParts.insert(PolyParts.end(), space.holes_begin(), space.holes_end());
            for (auto poly : PolyParts)
                for (int i = 0; i < poly.size(); i++)
                    pointIdVec.push_back({ poly[i],i });
            return pointIdVec;
        }

        std::string geojson;
        Context context;
        ControlProps controlProps;
    };
}

namespace SimplifyBoundary{
    SimplifyBoundaryTest::SimplifyBoundaryTest(const std::string& config_file){
        parseConfigFile(config_file, context, controlProps);
        geojson = parseInputFile(controlProps.InputFile);
        Solver ppSolver(geojson, context.simpProps.simplifyProps, context.simpProps.expandProps);
        showResult(ppSolver);
    }

    /// <summary>
    /// Parse intput file to string
    /// </summary>
    /// <param name="input_file_path"></param>
    /// <returns></returns>
    std::string SimplifyBoundaryTest::parseInputFile(const std::string& input_file_path) {
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
    void SimplifyBoundaryTest::parseConfigFile(const std::string& config_file, Context& context, ControlProps& controlProps) {
        Json::Value p = GeoJSON::parse_propsjson(config_file);
        std::string filename = p["filename"].asString();
        context.simpProps.simplifyProps.search_convex_area_thre = p["Sprops"]["search_convex_area_thre"].asDouble();
        context.simpProps.simplifyProps.search_concave_area_thre = p["Sprops"]["search_concave_area_thre"].asDouble();
        context.simpProps.simplifyProps.search_thre = p["Sprops"]["search_thre"].asDouble();
        context.simpProps.simplifyProps.del_area_thre = p["Sprops"]["del_area_thre"].asDouble();
        context.simpProps.simplifyProps.del_mdis_thre = p["Sprops"]["del_mdis_thre"].asDouble();
        context.simpProps.simplifyProps.del_bbs_width_thre = p["Sprops"]["del_bbs_width_thre"].asDouble();
        context.simpProps.simplifyProps.not_del_pbrate = p["Sprops"]["not_del_pbrate"].asDouble();
        context.simpProps.simplifyProps.theta_thre1 = p["Sprops"]["theta_thre1"].asDouble();
        context.simpProps.simplifyProps.theta_thre2 = p["Sprops"]["theta_thre2"].asDouble();
        context.simpProps.simplifyProps.min_circle_points_num = p["Sprops"]["min_circle_points_num"].asDouble();
        context.simpProps.simplifyProps.k1 = p["Sprops"]["k1"].asDouble();
        context.simpProps.simplifyProps.k2 = p["Sprops"]["k2"].asDouble();
        context.simpProps.simplifyProps.thre_max_multipe = p["Sprops"]["thre_max_multipe"].asDouble();
        context.simpProps.simplifyProps.thre_gap_node_num = p["Sprops"]["thre_gap_node_num"].asDouble();
        context.simpProps.simplifyProps.isRemerge = p["Sprops"]["isRemerge"].asBool();
        context.simpProps.simplifyProps.isProcess = p["Sprops"]["isProcess"].asBool();
        context.simpProps.expandProps.simplify_order = p["Eprops"]["simplify_order"].asString();
        context.simpProps.expandProps.offset = p["Eprops"]["offset"].asDouble();
        context.simpProps.expandProps.tri_simplify_cost = p["Eprops"]["tri_simplify_cost"].asDouble();
        context.simpProps.expandProps.isPostProcess = p["Eprops"]["isPostProcess"].asBool();
        controlProps.showId = p["Cprops"]["showId"].asBool();
        controlProps.showMode = p["Cprops"]["showMode"].asInt();
        controlProps.InputFile = p["Cprops"]["InputFile"].asString();
        controlProps.OutputFile = p["Cprops"]["OutputFile"].asString();
    }

    /// <summary>
    /// Show result with QT
    /// </summary>
    /// <param name="ppSolver"></param>
    void SimplifyBoundaryTest::showResult(const Solver& sbSolver) {
        std::cout << "Show Result:" << std::endl;
        const Polygon_with_holes_2& space = sbSolver.origin_space;
        std::vector<std::pair<Point_2, int>>& pointIdVec = get_point_id(space);
        const Polygon_with_holes_2& space_simplify = sbSolver.simplify_space;
        if (controlProps.OutputFile != "" && controlProps.OutputFile.substr(controlProps.OutputFile.length() - 7, 7) == "geojson") {
            std::string simpliy_geojson = GeoJSON::Pwh_to_geojson<K>(space_simplify);
            size_t lastSlashPos = controlProps.InputFile.find_last_of('/');
            if (lastSlashPos != std::string::npos) 
                write_geojson(simpliy_geojson, controlProps.OutputFile);
        }
        // cal bbs
        auto outer = space_simplify.outer_boundary();
        Vector_2 polygon_offset = Vector_2(outer.left_vertex()->x(), outer.bottom_vertex()->y());
        double polygon_width = CGAL::to_double(outer.right_vertex()->x() - outer.left_vertex()->x());
        double polygon_height = CGAL::to_double(outer.top_vertex()->y() - outer.bottom_vertex()->y());
        const std::vector<Point_2>& space_simplify_seg_point = sbSolver.seg_points;

        std::vector<Polygon_2> PolyParts_outer, PolyParts_holes;
        PolyParts_outer.push_back(space.outer_boundary());
        PolyParts_holes.insert(PolyParts_holes.end(), space.holes_begin(), space.holes_end());

        std::vector<Polygon_2> PolyParts_simplify_outer, PolyParts_simplify_holes;
        PolyParts_simplify_outer.push_back(space_simplify.outer_boundary());
        PolyParts_simplify_holes.insert(PolyParts_simplify_holes.end(), space_simplify.holes_begin(), space_simplify.holes_end());
        const std::vector<Polygon_2>& PolyParts_simplify_del_poly = sbSolver.del_polys;

        std::cout << "start gui" << std::endl;
        // # Qt gui show
        int argc = 1;
        const char* argv[2] = { "t2_viewer","\0" };
        QApplication app(argc, const_cast<char**>(argv));
        int screenWidth = 1920;
        int screenHeight = 1080;
        int windowWidth = 720;
        int windowHeight = 540;

        CGAL::Color holes_color(255, 255, 255);
        CGAL::Color segment_color(0, 0, 0);
        CGAL::Color face_color(67, 177, 235);
        CGAL::Color point_color(255, 0, 0);
        CGAL::Color del_segment_color(0, 0, 255);
        CGAL::Color del_face_color(0, 0, 50);
        CGAL::Color seg_point_color(255, 128, 0);
        CGAL::Color font_color(0, 255, 255);

        int centerX2 = screenWidth * 2 / 4;
        std::string title2 = "viewer-" + controlProps.InputFile;
        CGAL::CenterLineViewer<K> mainwindow2(app.activeWindow(), *PolyParts_simplify_outer.begin(), title2.c_str());
        if (controlProps.showId)
            for (auto p2i : pointIdVec) {
                mainwindow2.drawText(p2i.first, QString::number(p2i.second));
            }
        mainwindow2.drawPoints(space_simplify_seg_point, seg_point_color);
        mainwindow2.drawPartitions(PolyParts_simplify_del_poly, del_segment_color, del_face_color, point_color);
        mainwindow2.drawPartitions(PolyParts_simplify_holes, segment_color, holes_color, point_color);
        mainwindow2.drawPartitions(PolyParts_simplify_outer, segment_color, face_color, point_color);
        mainwindow2.move(centerX2 - mainwindow2.width() / 2, screenHeight / 2 - mainwindow2.height() / 2);

        std::cout << "PolyParts Size: outer(" << PolyParts_outer[0].size() << ")";
        if (PolyParts_holes.size() > 0)
            std::cout << " inner(" << PolyParts_holes[0].size() << ")";
        std::cout << endl;
        std::cout << "PolyParts_simply Size: outer(" << PolyParts_simplify_outer[0].size() << ")";
        if (PolyParts_simplify_holes.size() > 0)
            std::cout << " inner(" << PolyParts_simplify_holes[0].size() << ")";
        std::cout << std::endl;

        mainwindow2.show();
        app.exec();
    }

    /**
    * @brief Write a geojson to a file
    *
    * @param[in] geojson	A long string.
    * @param[in] outputPath Output Path.
    **/
    void SimplifyBoundaryTest::write_geojson(const std::string& geojson, const std::string& outputPath) {
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