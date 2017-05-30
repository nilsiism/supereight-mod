#ifndef TIMINGS_H
#define TIMINGS_H

void synchroniseDevices();

inline double tock() {
    synchroniseDevices();
#ifdef __APPLE__
		clock_serv_t cclock;
		mach_timespec_t clockData;
		host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
		clock_get_time(cclock, &clockData);
		mach_port_deallocate(mach_task_self(), cclock);
#else
    struct timespec clockData;
    clock_gettime(CLOCK_MONOTONIC, &clockData);
#endif
    return (double) clockData.tv_sec + clockData.tv_nsec / 1000000000.0;
}

#endif //TIMINGS_H
