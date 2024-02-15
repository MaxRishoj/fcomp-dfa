#ifndef LOG_H_
#define LOG_H_

// Module for logging information. We associate a category, so we can easily switch them on/off.

#include "common.h"

// NOTE: These are directly casted to ints, so do not assign them values manually.
const i32 kLogCount = 2;
enum Log { Stat, Debug };

void init_log();
void log_enable(const Log log, bool enabled = true);

void log(const Log log, const char *format, ...);

#endif // LOG_H_
