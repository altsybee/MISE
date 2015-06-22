#ifndef DISTANCEENTRY_H
#define DISTANCEENTRY_H


struct DistanceEntry
{
//        DistanceEntry& operator=(DistanceEntry& v);//перегрузка
    float dist;
    int x; // particle id from 1st array
    int y; // particle id from 2nd array
    bool inInteraction;
    friend bool operator > (DistanceEntry &c1, DistanceEntry &c2);
    friend bool operator < (DistanceEntry &c1, DistanceEntry &c2);
    friend bool operator >= (DistanceEntry &c1, DistanceEntry &c2);
    friend bool operator <= (DistanceEntry &c1, DistanceEntry &c2);
};


#endif // DISTANCEENTRY_H
