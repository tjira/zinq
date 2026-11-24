#pragma once

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/// Opens an interactive OpenGL viewport rendering the molecular system or trajectory at the given file path.
bool viewer_open_window(const char *filename);

#ifdef __cplusplus
}
#endif
