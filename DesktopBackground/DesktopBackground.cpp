// DesktopBackground.cpp : Defines the entry point for the application.
//

#define _CRT_SECURE_NO_WARNINGS
#include "models.h"


static constexpr std::uint32_t DEFAULT_WINDOWS_DPI = 96;
/* manifest tool cmd:
* cd C:\Users\mglas\source\repos\DesktopBackground\x64\Debug
* mt.exe -manifest ../../DesktopBackground/DesktopBackground.exe.manifest -outputresource:DesktopBackground.exe;1
*/

std::uint64_t get_current_time() {
     return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
}

template <typename T>
T random_int(T a, T b) {
    static std::random_device device{};
    static std::default_random_engine engine(device());
    const std::uniform_int_distribution<T> dist(a, b);
    return dist(engine);
}

template <typename T>
T random_real(T a, T b) {
    static std::random_device device{};
    static std::default_random_engine engine(device());
    std::uniform_real_distribution<T> dist(a, b);
    return dist(engine);
}

std::int32_t quadrant(double theta) {
    return (static_cast<std::int32_t>(theta * 2 / std::numbers::pi) % 4) + 1;
}


std::int32_t sign(int val) {
    if (val == 0) { return val; }
    return (val > 0) * 2 - 1;
}

/* this function copied and edited from https://gist.github.com/kingseva/a918ec66079a9475f19642ec31276a21 */
void BindStdHandlesToConsole() {
    //TODO: Add Error checking.
    
    // Redirect the CRT standard input, output, and error handles to the console
    freopen("CONIN$", "r", stdin);
    freopen("CONOUT$", "w", stderr);
    freopen("CONOUT$", "w", stdout);
    
    // Note that there is no CONERR$ file
    HANDLE hStdout = CreateFile(L"CONOUT$", GENERIC_READ|GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    HANDLE hStdin = CreateFile(L"CONIN$", GENERIC_READ|GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    
    SetStdHandle(STD_OUTPUT_HANDLE, hStdout);
    SetStdHandle(STD_ERROR_HANDLE, hStdout);
    SetStdHandle(STD_INPUT_HANDLE, hStdin);

    //Clear the error state for each of the C++ standard stream objects. 
    std::wclog.clear();
    std::clog.clear();
    std::wcout.clear();
    std::cout.clear();
    std::wcerr.clear();
    std::cerr.clear();
    std::wcin.clear();
    std::cin.clear();
}


/* This section (next two functions, EnumWindowsProc and get_wallpaper_window) copied from https://www.anycodings.com/1questions/2323944/drawing-on-the-desktop-background-as-wallpaper-replacement-windowsc */
BOOL CALLBACK EnumWindowsProc(HWND hwnd, LPARAM lParam) {
    HWND p = FindWindowEx(hwnd, NULL, L"SHELLDLL_DefView", NULL);
    HWND *ret = (HWND*)lParam;

    if (p) {
        // Gets the WorkerW Window after the current one.
        *ret = FindWindowEx(NULL, hwnd, L"WorkerW", NULL);
    }
    return true;
}

HWND get_wallpaper_window() {
    // Fetch the Progman window
    HWND progman = FindWindow(L"Progman", NULL);
    // Send 0x052C to  Progman. This message directs Progman to spawn a 
    // WorkerW behind the desktop icons. If it is already there, nothing 
    // happens.
    // Progman will only accept this if Animation Effects (or similar setting) is turned on
    SendMessageTimeout(progman, 0x052C, 0, 0, SMTO_NORMAL, 1000, nullptr);
    // We enumerate all Windows, until we find one, that has the SHELLDLL_DefView  
    // as a child. 
    // If  we found that window, we take its next sibling and assign it to workerw.
    HWND wallpaper_hwnd = nullptr;
    EnumWindows(EnumWindowsProc, (LPARAM)&wallpaper_hwnd);

    return wallpaper_hwnd;
}


int WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd) {
    AllocConsole();
    BindStdHandlesToConsole();
    std::ios_base::sync_with_stdio(false);
    sf::ContextSettings settings;
    settings.antialiasingLevel = 2;

    HWND wallpaper = get_wallpaper_window();
    RECT workarea;
    SystemParametersInfo(SPI_GETWORKAREA, NULL, &workarea, NULL);
    const int width = workarea.right - workarea.left, height = workarea.bottom - workarea.top;
    SetWindowPos(wallpaper, NULL, 0, 0, width, height, 0);
    
    sf::RenderWindow window(wallpaper, settings);
    if (!window.isOpen() || !sf::Shader::isAvailable()) {
        return 1;
    }

    /* shader files are copied and edited from https://en.sfml-dev.org/forums/index.php?topic=7827.0 */

    sf::Shader blooms;
    if (!blooms.loadFromFile("./bloom.fragment", sf::Shader::Fragment)) {
        std::cout << "bloom shader loading error\n";
    }

    const float sigma = 5.0F;
    const float glow_multiplier = 0.4F;
    
    blooms.setUniform("sourceTexture", sf::Shader::CurrentTexture);
    blooms.setUniform("sigma", sigma);
    blooms.setUniform("glowMultiplier", glow_multiplier);
    blooms.setUniform("width", static_cast<float>(width));
    blooms.setUniform("height", static_cast<float>(height));


#define PARTICLES

#ifdef PARTICLES

    const double display_scale = 1;
    const double ec = -10000;
    const double psize = 4;

    Simul s(width / (psize*2) + 1, height / (psize*2) + 1, 4); /* need extra cells in case dimensions are not nicely divisible */
    s.cellsize = psize*2;
    s.constraint_dim = {0.0, 0.0};
    s.constraint_sz = Vec2<double>(width, height);
    const sf::Vector2f constraint_point(s.constraint_dim.x * display_scale, s.constraint_dim.y * display_scale);

    const sf::Color bg = sf::Color(0xEA, 0xEA, 0xEA);
    const sf::Color fg = sf::Color(0x1C, 0x1C, 0x1C);
    const std::int32_t xoffset = (width - s.constraint_sz.x * display_scale) / 2;
    const std::int32_t yoffset = (height - s.constraint_sz.y * display_scale) / 2;

    const std::uint64_t start = get_current_time();
    std::uint64_t last_create = start;
    std::uint64_t pwait = 500'000'000ULL;


    /*
    for (std::int32_t i = 1; i < s.x-1; i += 4) {
        for (std::int32_t j = 1; j < s.y-1; j += 4) {
            for (std::int32_t k = 0; k < s.cellsize; k += psize * 2) {
                for (std::int32_t l = 0; l < s.cellsize; l += psize * 2) {
                    Particle *p = new Particle(Vec2<double>(random_real(i*s.cellsize + k - 1, i*s.cellsize + k + 1), j*s.cellsize + l));
                    p->size = psize;
                    p->temperature = random_real(0.0, 50.0);
                    s.cells[i][j].particles.push_back(p);
                }
            }
        }
    }
    */

    POINT cpoint;
    sf::RenderTexture thisf;
    sf::Sprite sprite;
    thisf.create(width, height, settings);
    sf::Font font;
    font.loadFromFile("./consola.ttf");
    sf::Text stext("", font);
    sf::Text rtext("", font);
    sf::Text otext("", font);
    rtext.setPosition(0, height - rtext.getCharacterSize());
    stext.setPosition(0, height - rtext.getCharacterSize() - stext.getCharacterSize());
    otext.setPosition(0, height - rtext.getCharacterSize() - stext.getCharacterSize() - otext.getCharacterSize());

    sf::Texture circletex;
    sf::Sprite point;
    circletex.loadFromFile("./circle.png");
    point.setTexture(circletex);
    const double initscale = circletex.getSize().x/2.0; /* we don't care x or y its circle */

    std::int32_t ocount = 0;
	std::uint64_t render_start = get_current_time();
    while (true) {
        /*
        std::uint64_t now = get_current_time();
        if (now - start < 30 * 1'000'000'000ULL) {
            if (now - last_create >= pwait) {
        */
        const std::uint64_t simul_start = get_current_time();
        if (ocount < 7000) {
            for (std::int32_t i = 0; i < 3; i++) {
                Particle *p = new Particle({random_real<double>(psize, width - psize), random_real<double>(psize, height/2.0 - psize)});
                p->resistance = random_real(0.12, 4.0);
                p->size = psize;
                p->temperature = random_real(0.0, 120.0);
                s.addparticle(p);
            }
        }
        s.update(1/60.0, 8);
        const double simul_length = (get_current_time() - simul_start) / 1'000'000.0L; /* in ms */

        /* std::this_thread::sleep_for(std::chrono::milliseconds(100/120)); */
        /* drawing cells */
        /*
        for (std::int32_t i = 0; i < s.x; i++) {
            for (std::int32_t j = 0; j < s.y; j++) {
                sf::RectangleShape cellr(sf::Vector2f(s.cellsize * display_scale, s.cellsize * display_scale));
                sf::RectangleShape cellbr(sf::Vector2f(s.cellsize * display_scale - 2, s.cellsize * display_scale - 2));
                cellr.setPosition(i * s.cellsize * display_scale + xoffset, j * s.cellsize * display_scale + yoffset);
                cellr.setFillColor(sf::Color::Blue);
                cellbr.setPosition(i * s.cellsize * display_scale + 1 + xoffset, j * s.cellsize * display_scale + 1 + yoffset);
                cellbr.setFillColor(fg);
                window.draw(cellr);
                window.draw(cellbr);
            }
        }
        */
        /*
        GetCursorPos(&cpoint);
        const Vec2<double> cposition(cpoint.x, cpoint.y);
        */
        ocount = 0;
        thisf.clear(fg);
        for (std::int32_t i = 0; i < s.x; i++) {
            for (std::int32_t j = 0; j < s.y; j++) {
                for (Particle *pparticle : s.cells[i][j].particles) {
                    Particle& particle = *pparticle;
                    ocount++;
                    if (particle.temperature < 80) { continue; }
                    point.setScale(sf::Vector2f(particle.size * display_scale / initscale, particle.size * display_scale / initscale));
                    const std::int32_t tempmul = std::clamp<std::int32_t>(particle.temperature, 0, 255);
                    point.setColor(sf::Color(tempmul, static_cast<std::int32_t>(std::pow(tempmul/255.0, 2)*255/3.0), static_cast<std::int32_t>(std::pow(tempmul/255.0, 3)*255/10.0))); 
                    point.setPosition((particle.pos_cur.x - particle.size) * display_scale + xoffset, (particle.pos_cur.y - particle.size) * display_scale + yoffset);
                    thisf.draw(point);
                    /*
                    const Vec2<double> cdiff(cposition - particle.pos_cur);
                    particle.accelerate((cdiff / cdiff.mod()) * ec * particle.size / std::clamp<double>(cdiff.mod()*cdiff.mod(), 0.1, 10000));
                    */
                }
            }
        }
        thisf.display();
        window.clear();
        sprite.setTexture(thisf.getTexture());
        window.draw(sprite, &blooms);
        //window.draw(sprite);

        stext.setString("simul: " + std::to_string(std::round(simul_length * 1000.0)/1000.0) + " ms");
        otext.setString("objects: " + std::to_string(ocount));
        window.draw(stext);
        window.draw(otext);

        const double render_length = (get_current_time() - render_start - simul_length) / 1'000'000.0L; /* in ms */
		render_start = get_current_time();
        rtext.setString("render: " + std::to_string(std::round(render_length * 1000.0)/1000.0) + " ms");
        window.draw(rtext);
        window.display();
    }

#else

    constexpr std::int32_t ring_padding = 5; /* px */

    std::random_device device{};
    std::default_random_engine engine(device());
    
    phyanim::Ring rings[5];
    
    const sf::Color colors[] = { sf::Color(127, 51, 40), sf::Color(197, 108, 75),
        sf::Color(219, 158, 44), sf::Color(129, 175, 160), sf::Color(92, 112, 126) };
    const sf::Color dimmer(30, 30, 30);

    /* init all rings */
    for (std::int32_t i = 0; i < 5; i++) {
        rings[i].width = std::uniform_int_distribution(5, 20)(engine);
        rings[i].top_color = colors[i];
        rings[i].bottom_color = sf::Color(colors[i].r - dimmer.r, colors[i].g - dimmer.g, colors[i].b - dimmer.b);
        rings[i].top_speed = std::uniform_real_distribution<double>(1, std::numbers::pi)(engine);
        rings[i].bottom_speed = rings[i].top_speed * 1.1;
        const std::int32_t this_top_arcs_count = std::uniform_int_distribution(3, 12)(engine);
        const std::int32_t this_bottom_arcs_count = this_top_arcs_count / 3;

        {
            rings[i].top.resize(this_top_arcs_count);
            std::vector<double> points(this_top_arcs_count * 2);
            for (std::int32_t j = 0 ; j < this_top_arcs_count * 2; j++) {
                points[j] = std::uniform_real_distribution<double>(0, 2 * std::numbers::pi)(engine);
            }
            std::sort(points.begin(), points.end());
            for (std::int32_t j = 0; j < this_top_arcs_count * 2; j += 2) {
                rings[i].top[j / 2].begin_angle = points[j];
                rings[i].top[j / 2].length = points[j + 1] - points[j];
            }
        }

        {
            rings[i].bottom.resize(this_bottom_arcs_count);
            std::vector<double> points(this_bottom_arcs_count * 2);
            for (std::int32_t j = 0 ; j < this_bottom_arcs_count * 2; j++) {
                points[j] = std::uniform_real_distribution<double>(0, 2 * std::numbers::pi)(engine);
            }
            std::sort(points.begin(), points.end());
            for (std::int32_t j = 0; j < this_bottom_arcs_count * 2; j += 2) {
                rings[i].bottom[j / 2].begin_angle = points[j];
                rings[i].bottom[j / 2].length = points[j + 1] - points[j];
            }
        }
    }
    
    sf::RenderTexture thisf;
    sf::Sprite sprite;
    thisf.create(width, height, settings);
    const std::int32_t points = 10;

    sf::VertexArray arcpoints(sf::Points);
	sf::VertexArray vertices(sf::TriangleFan, 2 * points);
    const std::int32_t xoffset = width / 2, yoffset = height / 2;
    while (true) {
        std::int32_t cur_dist = 100;
        for (std::int32_t ri = 0; ri < 5; ri++) {
            for (phyanim::Arc& arc : rings[ri].top) {
                arc.begin_angle += rings[ri].top_speed * 0.01;
                if (arc.begin_angle > 2 * std::numbers::pi) { arc.begin_angle -= 2 * std::numbers::pi; }

                /* TODO: ---- polar instead, this is too buggy ---- */

                double top_bound = std::tan(arc.begin_angle + arc.length);
                double bottom_bound = std::tan(arc.begin_angle);
                if (bottom_bound > top_bound) {
					std::cout << arc.begin_angle << ' ' << arc.length << std::endl;
                    double temp = top_bound;
                    top_bound = bottom_bound;
                    bottom_bound = temp;
                }
                /*
                if (quadrant(arc.begin_angle) != quadrant(arc.begin_angle + arc.length)) {
                    double temp = top_bound;
                    top_bound = bottom_bound;
                    bottom_bound = temp;
                }
                */

				const std::int32_t x_begin = cur_dist * std::cos(arc.begin_angle + arc.length);
				const std::int32_t y_begin = cur_dist * std::sin(arc.begin_angle);
                for (std::int32_t i = -cur_dist - rings[ri].width; i < (cur_dist + rings[ri].width); i++) {
                    for (std::int32_t j = -cur_dist - rings[ri].width; j < (cur_dist + rings[ri].width); j++) {
                        const double pdist = std::sqrt(i * i + j * j);
                        if (j < top_bound * i && j > bottom_bound * i && pdist > cur_dist && pdist < cur_dist + rings[ri].width) {
                            arcpoints.append(sf::Vertex(sf::Vector2f(i + xoffset, j + yoffset), rings[ri].top_color));
                        }
                    }
                }
            }
            cur_dist += ring_padding + rings[ri].width;
        }

        thisf.clear();
        thisf.draw(arcpoints);
        thisf.display();
        sprite.setTexture(thisf.getTexture());

        window.clear();
        window.draw(sprite);
        window.display();

        arcpoints.clear();
    }


#endif

    FreeConsole();

    return 0;
}
