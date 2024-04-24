#include "Partition.hpp"
#include "Placement.hpp"
#include "PartitionPlacement.hpp"
#include "PPGeoJSON.hpp"
#include "PartitionPlaceViewer.hpp"

#include <iostream>
#include <string>

namespace ParitionPlacement {
    struct PartitionPlacementTest
    {
        using PartitionProps = Partition::PartitionProps;
        using PlacementProps = Placement::PlacementProps;

        struct ControlProps {
            std::string InputFile;
            std::string OutputFile;
        };

        PartitionPlacementTest(const std::string& config_file);
        void parseConfigFile(const std::string& config_file, Context& context, ControlProps& controlProps);
        std::string parseInputFile(const std::string& input_file);
        void showResult(const Solver& ppSolver);

        std::string geojson;
    };

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
        Context context;
        ControlProps controlProps;
        parseConfigFile(config_file, context, controlProps);
        geojson = parseInputFile(controlProps.InputFile);
        Solver ppSolver(geojson, context);
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
    void PartitionPlacementTest::parseConfigFile(const std::string& config_file, Context& context, ControlProps& controlProps) {
        Json::Value p = parse_propsjson(config_file);
        controlProps.InputFile = p["ControlProps"]["InputFile"].asString();
        controlProps.OutputFile = p["ControlProps"]["OutputFile"].asString();
    }

    /// <summary>
    /// Show result with QT
    /// </summary>
    /// <param name="ppSolver"></param>
    void PartitionPlacementTest::showResult(const Solver& ppSolver) {
        std::cout << "Show Result:" << std::endl;
    }
}