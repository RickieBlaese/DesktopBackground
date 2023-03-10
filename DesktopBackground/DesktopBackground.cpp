// DesktopBackground.cpp : Defines the entry point for the application.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <thread>
#include <random>

#define NOMINMAX
#include "Windows.h"
#include "process.h"
#include "tlhelp32.h"
#include "shellscalingapi.h"
#include <SFML-2.5.1/include/SFML/Graphics.hpp>

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
    const std::uniform_real_distribution<T> dist(a, b);
    return dist(engine);
}

/* This section copied from https://www.anycodings.com/1questions/2323944/drawing-on-the-desktop-background-as-wallpaper-replacement-windowsc */
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


int main(int argc, char **argv) {
    AllocConsole();
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

    /* shader files are copied from https://en.sfml-dev.org/forums/index.php?topic=7827.0 */

    sf::Shader blooms;
    if (!blooms.loadFromFile("./bloom.fragment", sf::Shader::Fragment)) {
        std::cout << "bloom shader loading error\n";
    }

    const float sigma = 5.0F;
    const float glow_multiplier = 0.15F;
    
    blooms.setUniform("sourceTexture", sf::Shader::CurrentTexture);
    blooms.setUniform("sigma", sigma);
    blooms.setUniform("glowMultiplier", glow_multiplier);
    blooms.setUniform("width", static_cast<float>(width));
    blooms.setUniform("height", static_cast<float>(height));


    const double display_scale = 1;
    const double ec = -10000;
    const double psize = 4;

    Simul s(width / (psize*2) + 1, height / (psize*2) + 1); /* need extra in case dimensions are not nicely divisible */
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
    while (true) {
        /*
        std::uint64_t now = get_current_time();
        if (now - start < 30 * 1'000'000'000ULL) {
            if (now - last_create >= pwait) {
        */
        const std::uint64_t simul_start = get_current_time();
        if (ocount < 5000) {
            for (std::int32_t i = 0; i < 10; i++) {
                Particle *p = new Particle({random_real<double>(psize, width - psize), random_real<double>(psize, height/2.0 - psize)});
                p->size = psize;
                p->temperature = random_real(0.0, 50.0);
                s.addparticle(p);
            }
        }
        s.update(1/60.0, 4);
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
        const std::uint64_t render_start = get_current_time();
        GetCursorPos(&cpoint);
        const Vec2<double> cposition(cpoint.x, cpoint.y);
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
                    point.setColor(sf::Color(tempmul, tempmul/8, 0x05)); 
                    point.setPosition((particle.pos_cur.x - particle.size) * display_scale + xoffset, (particle.pos_cur.y - particle.size) * display_scale + yoffset);
                    thisf.draw(point);
                    const Vec2<double> cdiff(cposition - particle.pos_cur);
                    particle.accelerate((cdiff / cdiff.mod()) * ec * particle.size / std::clamp<double>(cdiff.mod()*cdiff.mod(), 0.1, 10000));
                }
            }
        }
        thisf.display();
        window.clear();
        sprite.setTexture(thisf.getTexture());
        window.draw(sprite, &blooms);

        stext.setString("simul: " + std::to_string(std::round(simul_length * 1000.0)/1000.0) + " ms");
        otext.setString("objects: " + std::to_string(ocount));
        window.draw(stext);
        window.draw(otext);

        const double render_length = (get_current_time() - render_start) / 1'000'000.0L; /* in ms */
        rtext.setString("render: " + std::to_string(std::round(render_length * 1000.0)/1000.0) + " ms");
        window.draw(rtext);
        window.display();
    }

    FreeConsole();

    return 0;
}
