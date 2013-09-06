// ConsoleColor.h


#pragma once
#include <iostream>
#ifdef WIN32
#include <windows.h>
#endif

inline std::ostream& blue(std::ostream &s)
{
#ifdef WIN32
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, FOREGROUND_BLUE
              |FOREGROUND_GREEN|FOREGROUND_INTENSITY);
#else
  system("Color 1");
#endif
    return s;
}

inline std::ostream& red(std::ostream &s)
{
#ifdef WIN32

    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, 
                FOREGROUND_RED|FOREGROUND_INTENSITY);
#else
  system("Color 4");
#endif

    return s;
}

inline std::ostream& green(std::ostream &s)
{
#ifdef WIN32

    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, 
              FOREGROUND_GREEN|FOREGROUND_INTENSITY);
#else
  system("Color 2");
#endif

    return s;
}

inline std::ostream& yellow(std::ostream &s)
{
#ifdef WIN32

    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, 
         FOREGROUND_GREEN|FOREGROUND_RED|FOREGROUND_INTENSITY);
#else
  system("Color 6");
#endif

    return s;
}

inline std::ostream& white(std::ostream &s)
{
#ifdef WIN32

    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout, 
       FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
#else
  system("Color 7");
#endif

    return s;
}
#ifdef WIN32
struct color {
    color(WORD attribute):m_color(attribute){};
    WORD m_color;
};

template <class _Elem, class _Traits>
std::basic_ostream<_Elem,_Traits>& 
      operator<<(std::basic_ostream<_Elem,_Traits>& i, color& c)
{

    HANDLE hStdout=GetStdHandle(STD_OUTPUT_HANDLE); 
    SetConsoleTextAttribute(hStdout,c.m_color);
    return i;
}
#endif
// Copyleft Vincent Godin
