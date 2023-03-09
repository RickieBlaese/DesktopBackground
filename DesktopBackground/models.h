#ifndef _MODELS_H
#define _MODELS_H

#include <vector>
#include <cmath>
#include <utility>
#include <cinttypes>


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

class Particle {
public:
    double temperature = 0, size = 1;
    Vec2<double> pos_old, pos_cur;
    Vec2<double> acceleration;


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
    /* x, y is number of cells, each 1x1 */
    std::int32_t x = 0, y = 0;
    double maxsize = 1;
    double cellsize = 1;
    Vec2<double> constraint_dim, constraint_sz; /* dim = x, y and sz = w, h */
    Cell **cells;

    Simul(std::int32_t x, std::int32_t y);
    ~Simul();

    void addparticle(Particle *pparticle);
    void accelerate(const Vec2<double>& acc);
    void check_particle(std::int32_t cx, std::int32_t cy, std::int32_t k);
    void update(double dt, std::int32_t substeps);

};


#endif
