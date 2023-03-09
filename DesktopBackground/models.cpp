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
	/* we don't actually care where since update will align it to the correct cell */
	cells[0][0].particles.push_back(pparticle);
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
	const double response_coeff = 0.75;
	const double elasticity = 0.2;
	accelerate({0.0, 15.0}); /* gravity */
	for (std::int32_t i = 0; i < x; i++) {
		for (std::int32_t j = 0; j < y; j++) {
			for (std::int32_t ip = 0; ip < cells[i][j].particles.size(); ip++) {
				std::vector<Particle*>& particles = cells[i][j].particles;
				Particle *poparticle = cells[i][j].particles[ip];
				Particle& oparticle = *poparticle;

				oparticle.update(dt);

				/* check surrounding */
				for (std::int32_t si = 0; si < substeps; si++) {
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
									const double delta = 0.5 * response_coeff * (dist - oparticle.size - cparticle.size);
									/* update positions */
									oparticle.pos_cur -= nv * (oratio * delta);
									cparticle.pos_cur += nv * (cratio * delta);
									oparticle.pos_cur += nv * 0.5 * elasticity;
									cparticle.pos_cur -= nv * 0.5 * elasticity;
									check_particle(n + i, m + j, ik);
								}
							}
						}
					}
				}
				/* circle bounds checking */
				const Vec2 v = constraint_center - oparticle.pos_cur;
				const double dist = oparticle.pos_cur.dist(constraint_center);
				if (dist > (crad - oparticle.size)) {
					const Vec2<double> nv = v / dist;
					oparticle.pos_cur = constraint_center - nv * (crad - oparticle.size);
				}
				/* have to re find in case a particle was moved out of cell due to earlier check_particle */
				check_particle(i, j, std::find(particles.begin(), particles.end(), poparticle) - particles.begin());
			}
		}
	}
}
