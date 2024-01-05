#ifndef _MODELS_H
#define _MODELS_H

#include <atomic>
#include <vector>
#include <utility>
#include <mutex>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <thread>
#include <numbers>
#include <random>

#include <cmath>
#include <cinttypes>

#define NOMINMAX
#include "Windows.h"
#include "process.h"
#include "tlhelp32.h"
#include "shellscalingapi.h"
#include <SFML-2.5.1/include/SFML/Graphics.hpp>

/* this function just divides it up then adds remainder to last of the sections */
void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out);


template <typename T>
struct Vec2 {
    T x = 0, y = 0;

    Vec2() = default;
    Vec2(const T& x, const T& y);

    Vec2 operator=(const Vec2& other) noexcept;
    Vec2 operator+(const Vec2& other) const noexcept;
    Vec2 operator+=(const Vec2& other) noexcept;
    Vec2 operator-() const noexcept;
    Vec2 operator-(const Vec2& other) const noexcept;
    Vec2 operator-=(const Vec2& other) noexcept;
    Vec2 operator*(const T& other) const noexcept;
    Vec2 operator/(const T& other) const;
    bool operator==(const Vec2<T>& other) const noexcept;

    double mod() const noexcept;
    double dist(const Vec2& other) const;
};

using MPos = Vec2<std::int32_t>;

class Lock {
public:
    std::atomic_int32_t i = 0;
    std::atomic_int32_t have_exited = 0;
    std::int32_t total = 0;

    Lock() = default;
    Lock(std::int32_t total);
    ~Lock() = default;

    void operator=(const Lock& other);

    void wait();
    void wait(std::stop_token stoken);

};

class Particle {
public:
    double temperature = 0, size = 1, resistance = 1;
    Vec2<double> pos_old, pos_cur;
    Vec2<double> acceleration;
    MPos cell;


    Particle() = default;
    Particle(Vec2<double> pos);
    ~Particle() = default;

    bool operator==(const Particle& other) const noexcept;

    void update(double dt);
    void accelerate(const Vec2<double>& acc);

};

class Cell {
public:
    MPos mpos;
    std::vector<Particle*> particles;

    Cell();
    Cell(MPos mpos);
    ~Cell();

};

class Simul {
public:
    /* x, y is number of cells */
    std::int32_t x = 0, y = 0;
    std::int32_t thread_count = 1;
    double maxsize = 1;
    double cellsize = 1;
    double temp_trans = 40; /* how much per second */
    double temp_decay = 90; /* per sec */
    double temp_wall_decay = 0;
    double temp_gain = 700;
    double elasticity = 0.75;
    Vec2<double> constraint_dim, constraint_sz; /* dim = x, y and sz = w, h */
    Cell **cells;
    std::vector<Particle*> particles;

    std::jthread *threads;
    /* section for each thread to mess with: bx, by are begin coordinates and x, y are end */
    std::atomic<double> dt = 0;
    std::vector<std::tuple<std::int32_t, std::int32_t, std::int32_t, std::int32_t>> tcoords; /* does not need to be atomic since lock controls */
    std::atomic_int32_t substeps = 0;
    std::mutex particle_mutex;
    Lock begin_lock;

    Simul(std::int32_t x, std::int32_t y, std::int32_t thread_count);
    ~Simul();

    void addparticle(Particle *pparticle);
    void accelerate(const Vec2<double>& acc);
    void check_particle(Particle *pparticle);
    void update(double dt, std::int32_t substeps);

};

namespace phyanim {
	struct Arc {
		double begin_angle = 0.0, length = 0.0;
	};

	struct Ring {
		std::int32_t width = 0.0;
		std::vector<Arc> top, bottom;
        double top_speed = 0.0, bottom_speed = 0.0;
		sf::Color top_color, bottom_color;
	};

} /* namespace phyanim */


#endif
