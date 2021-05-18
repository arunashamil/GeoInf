import requests
import math as m
import numpy as np
from pyorbital.orbital import Orbital
from datetime import datetime, timedelta
import matplotlib.pyplot as plt


def to_radians(value): #перевод в радианы
    return (value * np.pi) / 180


def to_degrees(value): #перевод в градусы
    return (value * 180) / np.pi


def to_decarts(r, theta, phi): #из сферических в декартовы координаты
    x = r * m.cos(theta) * m.cos(phi)
    y = r * m.cos(theta) * m.sin(phi)
    z = r * m.sin(theta)
    return x, y, z


def angle_V1_V2(V1, V2, V1_length, V2_length): #угол между векторами
    return m.acos((V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]) / (V1_length * V2_length))


def scalar_pr(V1, V2): #скалярное произведение векторов
    return V1[0] * V2[0] + V1[1] * V2[1] + V1[2] * V2[2]


def set_V(x1, y1, z1, x2, y2, z2): #задать вектор(его координаты и длину)
    V = [x2 - x1, y2 - y1, z2 - z1]
    V_length = m.sqrt(V[0] ** 2 + V[1] ** 2 + V[2] ** 2)
    return V, V_length


def get_satellite_data(tle_1, tle_2, utc_time): #получение данных со спутника в виде долготы, широты, высоты над Землёй
    orb = Orbital("N", line1=tle_1, line2=tle_2)
    lon, lat, height_st = orb.get_lonlatalt(utc_time)
    return lon, lat, height_st


def get_tle(file_with_tle, satellite_name): #получение TLE
    record = requests.get(file_with_tle, stream=True)
    open('TLE.txt', 'wb').write(record.text)
    file = open('TLE.txt', 'r')
    temporary = file.read().split("\n")[:-1]
    for i in range(len(temporary)):
        if temporary[i] == satellite_name:
            return [satellite_name, temporary[i + 1], temporary[i + 2]]


def dist_to_PL(x1, y1, z1, x2, y2, z2): #расстояние до касательной к сфере плоскости от точки вне сферы
    D = -x1 ** 2 - y1 ** 2 - z1 ** 2
    dif = (x1 * x2 + y1 * y2 + z1 * z2 + D) / ((x1 ** 2 + y1 ** 2 + z1 ** 2) ** 0.5)
    return dif


h_LK = 0.197 #высота ЛК над морем
R = 6378.1375 #радиус Земли
latitude_LK_R = to_radians(55.928895) #широта ЛК
longitude_LK_R = to_radians(37.521498) #долгота ЛК
x_LK, y_LK, z_LK = to_decarts(R + h_LK, latitude_LK_R,
                              longitude_LK_R) #координаты ЛК в декартовой системе
D = -(x_LK ** 2 + y_LK ** 2 + z_LK ** 2) #свободный член в уравнении касательной плоскости
z_P = -D / z_LK #z-ая координата пересечения плоскости с осью апликат Земли
y_P = -D / y_LK
x_P = -D / x_LK
North_V, North_V_L = set_V(0, 0, z_P, x_LK, y_LK, z_LK)  #вектор, указывающий на север в данной плоскости
Normal_V, Normal_V_L = set_V(0, 0, 0, x_LK, y_LK, z_LK) #вектор нормали к плоскости
East_V, East_V_L = set_V(0, 0, 0, North_V[1] * Normal_V[2] - Normal_V[1] * North_V[2], #вектор, указывающий на восток(через векторное произведение)
                         -(North_V[0] * Normal_V[2] - North_V[2] * Normal_V[0]),
                         North_V[0] * Normal_V[1] - Normal_V[0] * North_V[1])

#списки для отбора нужных данных:
array_x = []
array_y = []
array_z = []
array_time = []
array_elevation = []
array_azimuth = []
array_usual_elevation = []
array_usual_azimuth = []
array_usual_time = []
temporary_elevation_data = []
temporary_azimuth_data = []
temporary_time_data = []

list_tle = get_tle("https://celestrak.com/NORAD/elements/active.txt", "NOAA 19                 ")
print(list_tle)
print("Your date start MM HH DD MM YYYY: ")
#ввод начальной даты(local):
start_time = input().split()
start_time = [int(c) for c in start_time]
start_time = [int(c) for c in start_time]
start_time = datetime(start_time[4], start_time[3], start_time[2], start_time[1], start_time[0])
#ввод конечной даты(local):
print("Your date end MM HH DD MM YYYY: ")
end_time = input().split()
end_time = [int(c) for c in end_time]
end_time = datetime(end_time[4], end_time[3], end_time[2], end_time[1], end_time[0])
minutes = int((end_time - start_time).total_seconds() / 60) #разница между датами в минутах

for i in range(minutes):
    longitude, latitude, height = get_satellite_data(list_tle[1], list_tle[2], start_time - timedelta(hours=3)) #время тут переводим в utc
    start_time = start_time + timedelta(minutes=1)
    latitude_radians = to_radians(latitude)
    longitude_radians = to_radians(longitude)
    x, y, z = to_decarts(height + R, latitude_radians, longitude_radians)
    array_x.append(x)
    array_y.append(y)
    array_z.append(z)
    Sat_V, Sat_V_L = set_V(x_LK, y_LK, z_LK, x, y, z) #определение вектора, указывающего на спутник от ЛК и его длины
    dist_to_P = dist_to_PL(x_LK, y_LK, z_LK, x, y, z)
    dist_to_Sat = ((x - x_LK) ** 2 + (y - y_LK) ** 2 + (z - z_LK) ** 2) ** 0.5 #расстояние от ЛК до спутника
    elevation = m.asin(dist_to_P / dist_to_Sat) #тут может быть отрицательные значения, потому что в формуле расстояния точки от плоскости мы убрали модуль в числителе
    azimuth = 0
    if elevation >= 0:
        Sat_V_N = [(Normal_V[0] * scalar_pr(Normal_V, Sat_V)) / (Normal_V_L ** 2),
                   (Normal_V[1] * scalar_pr(Normal_V, Sat_V)) / (Normal_V_L ** 2),
                   (Normal_V[2] * scalar_pr(Normal_V, Sat_V)) / (Normal_V_L ** 2)] #проекция вектора спутника на вектор нормали
        Sat_V_Pr, Sat_V_Pr_L = set_V(0, 0, 0, -Sat_V[0] + Sat_V_N[0], -Sat_V[1] + Sat_V_N[1], -Sat_V[2] + Sat_V_N[2]) #проекция вектора спутника на плоскость
        azimuth = angle_V1_V2(Sat_V_Pr, North_V, North_V_L, Sat_V_Pr_L)
        if 3 * np.pi / 2 >= angle_V1_V2(East_V, Sat_V_Pr, East_V_L, Sat_V_Pr_L) > np.pi / 2:
            azimuth = 2 * np.pi - azimuth
    elevation = to_degrees(elevation)
    array_azimuth.append(azimuth)
    array_time.append(start_time + timedelta(hours=3)) #перевод в local
    array_elevation.append(elevation)

for i in range(len(array_time)):
    if 180 >= array_elevation[i] >= 0:
        temporary_elevation_data.append(array_elevation[i])
        temporary_azimuth_data.append(array_azimuth[i])
        temporary_time_data.append(array_time[i])
    else:
        if len(temporary_elevation_data) != 0:
            array_usual_time.append(temporary_time_data)
            array_usual_azimuth.append(temporary_azimuth_data)
            array_usual_elevation.append(temporary_elevation_data)
        temporary_elevation_data = [] #используем списки из списков,чтобы выводить данные из каждого списка в списки отдельно(иначе если выводить сплошняком, то всё будет соединено ломаными линиями)
        temporary_azimuth_data = []
        temporary_time_data = []

for i in range(len(array_usual_time)):
    print("[ Time: ", array_usual_time[i][0], "] [ Azimuth: ", to_degrees(array_usual_azimuth[i][0]), "] [ Elevation: ",
          max(array_usual_elevation[i]), "]")

print(array_usual_elevation)
print(array_usual_azimuth)

fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim(bottom=90, top=0)
for phi, theta in zip(array_usual_azimuth, array_usual_elevation):
    ax.plot(phi, theta)
fig.set_size_inches(7, 7)
plt.show()

sf = plt.figure()
ax = sf.add_subplot(111, projection='3d')
ax.plot(array_x, array_y, array_z)
ax.scatter(x_LK, y_LK, z_LK, color='black')
sf.set_size_inches(7, 7)
plt.show()
