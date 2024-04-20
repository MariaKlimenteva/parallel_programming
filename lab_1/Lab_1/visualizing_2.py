import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Считывание данных из .csv файла
data = pd.read_csv("solution_1.csv")

# Создание фигуры и 3D осей
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Построение 3D графика
x = data["x"]
y = data["t"]
z = data["u"]
ax.scatter(x, y, z, c=z, cmap="viridis")

# Настройка осей и заголовка
ax.set_xlabel("x")
ax.set_ylabel("t")
ax.set_zlabel("u(t, x)")
plt.title("Решение уравнения переноса (схема левый уголок)")

# Сохранение графика
plt.savefig("3d_graph_1.png")