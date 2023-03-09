// DesktopBackground.cpp : Defines the entry point for the application.
//

#include <iostream>
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

void kill_process(LPCWSTR filename) {
    HANDLE hSnapShot = CreateToolhelp32Snapshot(TH32CS_SNAPALL, NULL);
    PROCESSENTRY32 pEntry;
    pEntry.dwSize = sizeof(pEntry);
    BOOL hRes = Process32First(hSnapShot, &pEntry);
    while (hRes) {
        /* this is hacky but it works */
        if (strstr((char*)pEntry.szExeFile, (char*)filename) == 0) {
            HANDLE hProcess = OpenProcess(PROCESS_TERMINATE, 0,
                (DWORD)pEntry.th32ProcessID);
            if (hProcess != NULL) {
                TerminateProcess(hProcess, 9);
                CloseHandle(hProcess);
            }
        }
        hRes = Process32Next(hSnapShot, &pEntry);
    }
    CloseHandle(hSnapShot);
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

/* copied and edited from https://stackoverflow.com/questions/38334081/how-to-draw-circles-arcs-and-vector-graphics-in-sdl */
// rounding helper, simplified version of the function I use
int roundUpToMultipleOfEight(int v) {
    return (v + (8 - 1)) & -8;
}


int main(int argc, char **argv) {
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;

    HWND wallpaper = get_wallpaper_window();
    const int width = GetSystemMetrics(SM_CXVIRTUALSCREEN), height = GetSystemMetrics(SM_CYVIRTUALSCREEN);
    SetWindowPos(wallpaper, NULL, 0, 0, width, height, 0);
    
    sf::RenderWindow window(wallpaper);
    if (!window.isOpen()) {
        return 1;
    }



    const double display_scale = 3;

    Simul s(25, 25);
    s.cellsize = 8;
    s.constraint_center = {100.0, 100.0};
    s.crad = 100;
    const double psize = 2;
    const sf::Vector2f constraint_point(s.constraint_center.x, s.constraint_center.y);

    const sf::Color bg = sf::Color::Color(0xEA, 0xEA, 0xEA);
    const sf::Color fg = sf::Color::Color(0x1C, 0x1C, 0x1C);
    const std::int32_t xoffset = (width - 200 * display_scale) / 2;
    const std::int32_t yoffset = (height - 200 * display_scale) / 2;

    const std::uint64_t start = get_current_time();
    std::uint64_t last_create = start;
    std::uint64_t pwait = 500'000'000ULL;
    while (true) {
        std::uint64_t now = get_current_time();
        if (now - start < 10 * 1'000'000'000ULL) {
            if (now - last_create >= pwait) {
				Particle *p = new Particle({50, 130});
				p->size = psize;
				s.addparticle(p);
                last_create = now;
                pwait = random_int(5'000'000ULL, 10'000'000ULL);
            }
        }
		s.update(1/30.0);
        std::this_thread::sleep_for(std::chrono::milliseconds(100/120));
        window.clear(bg);
		sf::CircleShape bound(s.crad * display_scale);
        bound.setPosition(xoffset, yoffset);
        bound.setPointCount(100);
        bound.setFillColor(fg);
        window.draw(bound);
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
		for (std::int32_t i = 0; i < s.x; i++) {
			for (std::int32_t j = 0; j < s.y; j++) {
                for (const Particle *pparticle : s.cells[i][j].particles) {
                    const Particle& particle = *pparticle;
                    sf::CircleShape point(particle.size * display_scale);
                    point.setFillColor(sf::Color::Red);
                    point.setPosition((particle.pos_cur.x - particle.size) * display_scale + xoffset, (particle.pos_cur.y - particle.size) * display_scale + yoffset);
                    window.draw(point);
                }
            }
        }
        window.display();
    }


    return 0;
}
