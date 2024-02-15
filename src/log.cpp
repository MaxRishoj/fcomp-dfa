#include "log.h"

#include <cassert>
#include <cstdio>
#include <stdarg.h>

bool log_enabled[kLogCount];

void init_log() {
    // We disable all logs initially.
    for (i32 i = 0; i < kLogCount; i++) {
        log_enabled[i] = false;
    }
}

void log_enable(const Log log, bool enabled) {
    log_enabled[(i32)log] = enabled;
}

void log(const Log log, const char *format, ...) {
    if (!log_enabled[(i32)log]) return;

    switch (log) {
        case Log::Stat:  printf("[Stat]  "); break;
        case Log::Debug: printf("[Debug] "); break;
        default: assert(false);
    }

    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}
