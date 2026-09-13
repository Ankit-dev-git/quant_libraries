#include "DataContainerCreation.h"
#include <fstream>
#include <string>

void savedatacontainer(Grid& dc, const std::string& filename)
{
    std::ofstream file(filename);

    std::string buffer;
    buffer.reserve(1024 * 1024); // 1 MB

    for (int i = 0; i < dc.rows; ++i) {
        for (int j = 0; j < dc.cols; ++j) {

            buffer += std::to_string(dc(i, j));

            if (j < dc.cols - 1)
                buffer += ',';
        }

        buffer += '\n';

        if (buffer.size() >= 1024 * 1024) {
            file.write(buffer.data(), buffer.size());
            buffer.clear();
        }
    }

    // Write remaining data
    if (!buffer.empty())
        file.write(buffer.data(), buffer.size());
}