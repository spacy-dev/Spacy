#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace Spacy
{
    std::vector<double> lookupTableToVector(const std::string& line)
    {
        auto iss = std::istringstream(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>(iss),
                                        std::istream_iterator<std::string>()};
        std::vector<double> values;
        std::transform(begin(tokens), end(tokens),
                       std::back_inserter(values),
                       [](const std::string& token)
        {
            double value = std::stod(token);
            return value;
        });
        return values;
    }

    std::vector<double> extractLookupTable(const std::string& filename)
    {
        auto file = std::ifstream(filename);
        auto foundLookupTable = false;
        while(!file.eof())
        {
            std::string line;
            std::getline(file, line);
            if(foundLookupTable)
                return lookupTableToVector(line);
            foundLookupTable = line.compare(0, 12, "LOOKUP_TABLE") == 0;
        }
        return {};
    }

    bool compareLookupTable(const std::string& lhsVTKFileName, const std::string& rhsVTKFileName)
    {
        const auto lhsValues = extractLookupTable(lhsVTKFileName);
        const auto rhsValues = extractLookupTable(rhsVTKFileName);

        if(lhsValues.size() != rhsValues.size())
            return false;

        for(auto i=0u; i<lhsValues.size(); ++i)
            if(std::abs(lhsValues[i] - rhsValues[i]) > std::numeric_limits<double>::epsilon())
                return false;
        return true;
    }
}
