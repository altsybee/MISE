#include "DistanceEntry.h"

//DistanceEntry::DistanceEntry()
//{
//}

bool operator > (DistanceEntry &c1, DistanceEntry &c2)
{
    return c1.dist > c2.dist;
}
bool operator < (DistanceEntry &c1, DistanceEntry &c2)
{
    return c1.dist < c2.dist;
}
bool operator >= (DistanceEntry &c1, DistanceEntry &c2)
{
    return c1.dist >= c2.dist;
}
bool operator <= (DistanceEntry &c1, DistanceEntry &c2)
{
    return c1.dist <= c2.dist;
}
//DistanceEntry& DistanceEntry::operator=(DistanceEntry& v)//перегрузка
//{
//  x=v.x;y=v.y;dist=v.dist;
//  return *this;//возвращаем ссылку на текущий объект
//}
