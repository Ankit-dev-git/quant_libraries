#pragma once
#include <functional>

void parallel_for(int start, int end, const std::function<void(int, int)>& worker);
