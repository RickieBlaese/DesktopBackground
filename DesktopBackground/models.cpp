#include "models.h"


void partition(std::int32_t a, std::int32_t b, std::int32_t count, std::vector<std::pair<std::int32_t, std::int32_t>>& out) {
    if (count <= 1) {
        out.emplace_back(a, b);
        return;
    }
    const std::int32_t per = (b - a) / count;

    std::int32_t i = 0;
    for (; i < count - 1; i++) {
        out.emplace_back(i * per + a, (i + 1) * per + a);
    }
    out.emplace_back(i * per + a, b);
}

Lock::Lock(std::int32_t total) : total(total) {}

void Lock::operator=(const Lock& other) {
    i = other.i.load();
    total = other.total;
}

void Lock::wait() {
    i++;
    if (i >= total) {
        return;
    }
    /* this theoretically would have better branch prediction than the while (i < total) */
    while (true) {
        if (i >= total) { break; }
    }
    have_exited++;
    if (have_exited >= total) {
        i = 0;
        have_exited = 0;
    }
}

void Lock::wait(std::stop_token stoken) {
    i++;
    if (i >= total) {
        return;
    }
    /* this theoretically would have better branch prediction than the while (i < total) */
    while (!stoken.stop_requested()) {
        if (i >= total) { break; }
    }
    have_exited++;
    if (have_exited >= total) {
        i = 0;
        have_exited = 0;
    }
}



template <typename T>
Vec2<T>::Vec2(const T& x, const T& y) : x(x), y(y) {}

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



Simul::Simul(std::int32_t x, std::int32_t y, std::int32_t thread_count) : x(x), y(y), thread_count(thread_count) {
    cells = new Cell*[x];
    for (std::int32_t i = 0; i < x; i++) {
        cells[i] = new Cell[y];
        for (std::int32_t j = 0; j < y; j++) {
            cells[i][j].mpos = {i, j};
        }
    }

    begin_lock = Lock(thread_count + 1);
    std::vector<std::pair<std::int32_t, std::int32_t>> out;
    partition(0, x, thread_count * 2, out);
    for (std::int32_t i = 0; i < thread_count * 2; i++) {
        tcoords.emplace_back(out[i].first, 0, out[i].second, y);
    }

    auto update_mt = [&](std::stop_token stoken, std::int32_t index) {
        std::int32_t& fbx = std::get<0>(tcoords[index]);
        std::int32_t& fby = std::get<1>(tcoords[index]);
        std::int32_t& fx  = std::get<2>(tcoords[index]);
        std::int32_t& fy  = std::get<3>(tcoords[index]);

        std::int32_t& obx = std::get<0>(tcoords[index + 1LL]);
        std::int32_t& oby = std::get<1>(tcoords[index + 1LL]);
        std::int32_t& ox  = std::get<2>(tcoords[index + 1LL]);
        std::int32_t& oy  = std::get<3>(tcoords[index + 1LL]);
        const std::tuple<std::int32_t&, std::int32_t&, std::int32_t&, std::int32_t&> tcoord[] = {{fbx, fby, fx, fy}, {obx, oby, ox, oy}};
        
        while (!stoken.stop_requested()) {
            begin_lock.wait(stoken);
            for (auto& [bx, by, x, y] : tcoord) {
                for (std::int32_t i = bx; i < x; i++) {
                    for (std::int32_t j = by; j < y; j++) {
                        std::lock_guard guard(particle_mutex);
                        std::vector<Particle*>& particles = cells[i][j].particles;
                        for (std::int32_t ip = 0; ip < cells[i][j].particles.size(); ip++) {
                            Particle *poparticle = cells[i][j].particles[ip];
                            Particle& oparticle = *poparticle;
                            oparticle.temperature -= temp_decay * dt / substeps;
                            oparticle.accelerate(Vec2<double>(0, -oparticle.temperature * 40.0 / substeps));
                            oparticle.update(dt / substeps);

                            /* check surrounding */
                            for (std::int32_t n = -1; n < 2; n++) {
                                for (std::int32_t m = -1; m < 2; m++) {
                                    if (n + i < 0 || n + i >= x || m + j < 0 || m + j >= y) { continue; } /* out of bounds */
                                    for (std::int32_t ik = 0; ik < cells[n + i][m + j].particles.size(); ik++) {
                                        Particle *pcparticle = cells[n + i][m + j].particles[ik];
                                        if ((n == 0 && m == 0) && pcparticle == poparticle) { continue; } /* is itself */
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
                                            oparticle.temperature += (midtemp - oparticle.temperature) * temp_trans * dt / substeps / oparticle.resistance;
                                            cparticle.temperature += (midtemp - cparticle.temperature) * temp_trans * dt / substeps / oparticle.resistance;
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
                                oparticle.temperature += temp_gain * dt / substeps / oparticle.resistance;
                            }
                            oparticle.temperature = std::clamp(oparticle.temperature, 0.0, std::numeric_limits<double>::max());
                        }
                    }
                }
            }
        }
    };
    threads = new std::jthread[thread_count];
    for (std::int32_t i = 0; i < thread_count; i++) {
        threads[i] = std::jthread(update_mt, i * 2);
    }
}

Simul::~Simul() {
    for (std::int32_t i = 0; i < x; i++) {
        delete[] cells[i];
    }
    delete[] cells;
    delete[] threads; /* jthread does the signaling and joining for us in destructor */
    for (auto pparticle : particles) {
        delete pparticle;
    }
}

void Simul::addparticle(Particle *pparticle) {
    std::lock_guard guard(particle_mutex);
    particles.push_back(pparticle);
	std::int32_t mx = std::clamp(static_cast<std::int32_t>(pparticle->pos_cur.x) / static_cast<std::int32_t>(cellsize), 0, x - 1);
	std::int32_t my = std::clamp(static_cast<std::int32_t>(pparticle->pos_cur.y) / static_cast<std::int32_t>(cellsize), 0, y - 1);
	pparticle->cell = {mx, my};
    Cell &c = cells[mx][my];
	c.particles.push_back(pparticle);
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

void Simul::check_particle(Particle *pparticle) {
    Particle& particle = *pparticle;
    std::int32_t cx = particle.cell.x;
    std::int32_t cy = particle.cell.y;
    std::vector<Particle*>& particles = cells[cx][cy].particles;

    std::int32_t newx = particle.pos_cur.x / cellsize, newy = particle.pos_cur.y / cellsize;
    if (newx != cx || newy != cy) {
        particle_mutex.lock();
        auto ploc = std::find(particles.begin(), particles.end(), pparticle);
        if (ploc == particles.end()) {
            particle_mutex.unlock();
			return;
        }
        if (newx < 0 || newx >= x || newy < 0 || newy >= y) { /* out of bounds */
            particles.erase(ploc);
            this->particles.erase(std::find(this->particles.begin(), this->particles.end(), pparticle));
            particle_mutex.unlock();
            return;
        }
        cells[newx][newy].particles.push_back(pparticle);
        
        particles.erase(ploc);
        pparticle->cell = {newx, newy};
        particle_mutex.unlock();
    }
}

void Simul::update(double dt, std::int32_t substeps) {
    this->dt = dt;
    this->substeps = substeps;
    accelerate({0.0, 500.0}); /* gravity */
    for (std::int32_t si = 0; si < substeps; si++) {
        for (std::int32_t i = 0; i < thread_count; i++) {
            begin_lock.wait(); /* this will cause all the threads to stop waiting and then to update, since we initialized the lock to thread_count + 1 as the total */
        }
        std::lock_guard guard(particle_mutex);
        for (std::uint32_t i = 0; i < x; i++) {
            for (std::uint32_t j = 0; j < y; j++) {
                cells[i][j].particles.clear();
            }
        }
        for (std::uint32_t i = 0; i < particles.size(); i++) {
            Particle *pparticle = particles[i];
            std::int32_t mx = std::clamp(static_cast<std::int32_t>(pparticle->pos_cur.x) / static_cast<std::int32_t>(cellsize), 0, x - 1);
            std::int32_t my = std::clamp(static_cast<std::int32_t>(pparticle->pos_cur.y) / static_cast<std::int32_t>(cellsize), 0, y - 1);
            pparticle->cell = {mx, my};
            Cell &c = cells[mx][my];
            c.particles.push_back(pparticle);
            []{}();
        }
    }
}
