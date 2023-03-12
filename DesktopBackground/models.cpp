#include "models.h"


template <typename T>
Vec2<T>::Vec2<T>(const T& x, const T& y) : x(x), y(y) {}

template <typename T>
Vec2<T> Vec2<T>::operator=(const Vec2<T>& other) noexcept {
    x = other.x;
    y = other.y;
    return *this;
}

template <typename T>
Vec2<T> Vec2<T>::operator+(const Vec2<T>& other) const noexcept {
    return {x + other.x, y + other.y};
}

template <typename T>
Vec2<T> Vec2<T>::operator+=(const Vec2<T>& other) noexcept {
    return (*this = *this + other);
}

template <typename T>
Vec2<T> Vec2<T>::operator-() const noexcept {
    return {-x, -y};
}

template <typename T>
Vec2<T> Vec2<T>::operator-(const Vec2<T>& other) const noexcept {
    return *this + (-other);
}

template <typename T>
Vec2<T> Vec2<T>::operator-=(const Vec2<T>& other) noexcept {
    return (*this = *this + (-other));
}

template <typename T>
double Vec2<T>::mod() const noexcept {
    return std::sqrt(static_cast<double>(x * x) + static_cast<double>(y * y));
}

template <typename T>
double Vec2<T>::dist(const Vec2<T>& other) const {
    return (*this - other).mod();
}

template <typename T>
Vec2<T> Vec2<T>::operator*(const T& other) const noexcept {
    return {x * other, y * other};
}

template <typename T>
Vec2<T> Vec2<T>::operator/(const T& other) const {
    return {x / other, y / other};
}

template <typename T>
bool Vec2<T>::operator==(const Vec2<T>& other) const noexcept {
    return x == other.x && y == other.y;
}



Particle::Particle(Vec2<double> pos) : pos_old(std::move(pos)), pos_cur(pos_old) {}

bool Particle::operator==(const Particle& other) const noexcept {
    return size == other.size && temperature == other.temperature &&
        pos_old == other.pos_old && pos_cur == other.pos_cur && acceleration == other.acceleration;
}

void Particle::update(double dt) {
    const Vec2<double> disp = pos_cur - pos_old;
    pos_old = pos_cur;
    pos_cur = pos_cur + disp + acceleration * (dt * dt);
    acceleration = {0, 0};
}

void Particle::accelerate(const Vec2<double>& acc) {
    acceleration += acc;
}


Cell::Cell() : mpos(0, 0) {}

Cell::Cell(MPos mpos) : mpos(std::move(mpos)) {}

Cell::~Cell() {
    for (Particle *pparticle : particles) {
        delete pparticle;
    }
}



Simul::Simul(std::int32_t x, std::int32_t y) : x(x), y(y) {
    cells = new Cell*[x];
    for (std::int32_t i = 0; i < x; i++) {
        cells[i] = new Cell[y];
        for (std::int32_t j = 0; j < y; j++) {
            cells[i][j].mpos = {i, j};
        }
    }
}

Simul::~Simul() {
    for (std::int32_t i = 0; i < x; i++) {
        delete[] cells[i];
    }
    delete[] cells;
}

void Simul::addparticle(Particle *pparticle) {
    cells[1][1].particles.push_back(pparticle);
    check_particle(1, 1, cells[1][1].particles.size() - 1);
    maxsize = std::max(maxsize, pparticle->size);
}

void Simul::accelerate(const Vec2<double>& acc) {
    for (std::int32_t i = 0; i < x; i++) {
        for (std::int32_t j = 0; j < y; j++) {
            for (Particle *pparticle : cells[i][j].particles) {
                pparticle->accelerate(acc);
            }
        }
    }
}

void Simul::check_particle(std::int32_t cx, std::int32_t cy, std::int32_t k) {
    std::vector<Particle*>& particles = cells[cx][cy].particles;
    Particle& particle = *particles[k];

    std::int32_t newx = particle.pos_cur.x / cellsize, newy = particle.pos_cur.y / cellsize;
    if (newx != cells[cx][cy].mpos.x || newy != cells[cx][cy].mpos.y) {
        if (newx < 0 || newx >= x || newy < 0 || newy >= y) { /* out of bounds */
            particles.erase(particles.begin() + k);
            k--;
            return;
        }
        cells[newx][newy].particles.push_back(particles[k]);
        particles.erase(particles.begin() + k);
        k--;
    }
}

void Simul::update(double dt, std::int32_t substeps) {
    const double temp_trans = 20; /* how much per second */
    const double temp_decay = 80; /* per sec */
    const double temp_wall_decay = 100;
    const double temp_gain = 500;
    const double elasticity = 0.75;
    accelerate({0.0, 500.0}); /* gravity */
    for (std::int32_t si = 0; si < substeps; si++) {
        for (std::int32_t i = 0; i < x; i++) {
            for (std::int32_t j = 0; j < y; j++) {
                std::vector<Particle*>& particles = cells[i][j].particles;
                for (std::int32_t ip = 0; ip < cells[i][j].particles.size(); ip++) {
                    Particle *poparticle = cells[i][j].particles[ip];
                    Particle& oparticle = *poparticle;
                    oparticle.temperature -= temp_decay * dt / substeps;
                    oparticle.accelerate(Vec2<double>(0, -oparticle.temperature * 10.0 / substeps));
                    oparticle.update(dt / substeps);

                    /* check surrounding */
                    for (std::int32_t n = -1; n < 2; n++) {
                        for (std::int32_t m = -1; m < 2; m++) {
                            if (n + i < 0 || n + i >= x || m + j < 0 || m + j >= y) { continue; } /* out of bounds */
                            for (std::int32_t ik = 0; ik < cells[n + i][m + j].particles.size(); ik++) {
                                Particle *pcparticle = cells[n + i][m + j].particles[ik];
                                if ((n == 0 && m ==0) && pcparticle == poparticle) { continue; } /* is itself */
                                Particle& cparticle = *pcparticle;
                                if (cparticle.pos_cur.dist(oparticle.pos_cur) < oparticle.size + cparticle.size) {
                                    const Vec2<double> v = oparticle.pos_cur - cparticle.pos_cur;
                                    const double dist = oparticle.pos_cur.dist(cparticle.pos_cur);
                                    const Vec2<double> nv = v / dist;
                                    const double oratio = oparticle.size / (oparticle.size + cparticle.size);
                                    const double cratio = cparticle.size / (oparticle.size + cparticle.size);
                                    const double delta = 0.5 * (dist - (oparticle.size + cparticle.size));
                                    /* update positions */
                                    oparticle.pos_cur -= nv * (oratio * delta) * elasticity;
                                    cparticle.pos_cur += nv * (cratio * delta) * elasticity;
                                    const double midtemp = (oparticle.temperature + cparticle.temperature) / 2.0;
                                    oparticle.temperature += (midtemp - oparticle.temperature) * temp_trans * dt / substeps;
                                    cparticle.temperature += (midtemp - cparticle.temperature) * temp_trans * dt / substeps;
                                    check_particle(n + i, m + j, ik); /* cparticle */
                                }
                            }
                        }
                    }
                    /* circle bounds checking */
                    /*
                    const Vec2 v = constraint_center - oparticle.pos_cur;
                    const double dist = oparticle.pos_cur.dist(constraint_center);
                    if (dist > (crad - oparticle.size)) {
                        oparticle.temperature += temp_gain * dt;
                        const Vec2<double> nv = v / dist;
                        oparticle.pos_cur = constraint_center - nv * (crad - oparticle.size);
                    }
                    */
                    if (oparticle.pos_cur.x > constraint_dim.x + constraint_sz.x - oparticle.size) {
                        const Vec2<double> v(0, constraint_sz.y/2.0);
                        oparticle.pos_cur.x = constraint_dim.x + constraint_sz.x - oparticle.size;
                        oparticle.temperature -= temp_wall_decay * dt / substeps;
                    } else if (oparticle.pos_cur.x < constraint_dim.x + oparticle.size) {
                        oparticle.pos_cur.x = constraint_dim.x + oparticle.size;
                        oparticle.temperature -= temp_wall_decay * dt / substeps;
                    } else if (oparticle.pos_cur.y > constraint_dim.y + constraint_sz.y - oparticle.size) { /* is bottom */
                        oparticle.pos_cur.y = constraint_dim.y + constraint_sz.y - oparticle.size;
                        /*
                        if (oparticle.pos_cur.x > constraint_dim.x + constraint_sz.x/4.0 && oparticle.pos_cur.x < constraint_dim.x + constraint_sz.x * 3.0/4.0) {
                        }
                        */
                    } else if (oparticle.pos_cur.y < constraint_dim.y + oparticle.size) {
                        oparticle.pos_cur.y = constraint_dim.y + oparticle.size;
                        oparticle.temperature -= temp_wall_decay * dt / substeps;
                    } else if (oparticle.pos_cur.y > constraint_dim.y + constraint_sz.y * 0.98 - oparticle.size) {
                        oparticle.temperature += temp_gain * dt / substeps;
                    }
                    /* have to re find in case a particle was moved out of cell due to earlier check_particle */
                    check_particle(i, j, std::find(particles.begin(), particles.end(), poparticle) - particles.begin());
                }
            }
        }
    }
}
