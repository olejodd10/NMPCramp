#include "Timer.h"

#include <time.h>

#define TIMER_CLOCK CLOCK_REALTIME

static struct timespec then, now;

void timer_reset(void) {
    clock_gettime(TIMER_CLOCK, &then);
}

long timer_elapsed_ns(void) {
    clock_gettime(TIMER_CLOCK, &now);
    return (now.tv_sec-then.tv_sec)*1000000000 + now.tv_nsec-then.tv_nsec;
}
