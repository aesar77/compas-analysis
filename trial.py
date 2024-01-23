import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# a = np.array(0,dtype='bool')
# b = np.array([True, False, False, True], dtype='bool')

# print(np.shape(a), np.shape(b))
# c = np.concatenate([a,b])

# print(c)

# d = (np.concatenate([c,b]))
# print(d)
matplotlib.rcParams['figure.figsize'] = (15,10)
matplotlib.rcParams['lines.markersize'] = 1
matplotlib.rcParams['font.size'] = 14
rng = np.random.default_rng(seed=123) # random number generator
# size = 10000
# x1 = rng.normal(loc=5, scale=10, size=size)
# y1 = rng.normal(loc=5, scale=2, size=size)
# x2 = rng.normal(loc=-30, scale=5, size=size)
# y2 = rng.normal(loc=-20, scale=5, size=size)

# # Define hexbin grid extent
# xmin = min(*x1, *x2)
# xmax = max(*x1, *x2)
# ymin = min(*y1, *y2)
# ymax = max(*y1, *y2)
# xbin = np.linspace(xmin, xmax, 50)
# ybin = np.linspace(ymin, ymax, 50)

# # Draw figure with colorbars
# plt.figure(figsize=(10, 6))
# plt.hist2d(x1, y1, bins=(xbin, ybin) ,cmap='Blues')
# plt.colorbar(orientation='vertical')
# plt.hist2d(x2, y2,bins=(xbin, ybin), cmap='Reds', vmin=20)
# plt.colorbar(orientation='vertical')
# plt.show()

# hist1 = plt.hist2d(x1, y1, bins=(xbin, ybin) ,cmap='Blues')
# plt.colorbar(orientation='vertical')
# plt.show()
# hist2 = plt.hist2d(x2, y2,bins=(xbin, ybin), cmap='Reds')
# plt.colorbar(orientation='vertical')
# plt.show()

rng = np.random.default_rng(seed=123) # random number generator
size = 10000
x1 = rng.normal(loc=40, scale=5, size=size)
y1 = rng.normal(loc=40, scale=5, size=size)
x2 = rng.normal(loc=60, scale=5, size=size)
y2 = rng.normal(loc=60, scale=5, size=size)

# Define hexbin grid extent
xmin = min(*x1, *x2)
xmax = max(*x1, *x2)
ymin = min(*y1, *y2)
ymax = max(*y1, *y2)
ext = (xmin, xmax, ymin, ymax)

# Draw figure with colorbars
# plt.subplot(3,3,1)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.yscale("log")
# plt.xscale("log")
# plt.colorbar(hist1, orientation='vertical' )
# plt.colorbar(hist2, orientation='vertical') 

# plt.subplot(3,3,2)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.yscale("log")
# plt.xscale("log")
# plt.colorbar(hist1, orientation='vertical' )
# plt.colorbar(hist2, orientation='vertical')

# plt.subplot(3,3,3)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.colorbar(hist1, orientation='vertical', label="all" )
# plt.colorbar(hist2, orientation='vertical', label="undergoing")
# plt.yscale("log")
# plt.xscale("log") 

# plt.subplot(3,3,4)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.yscale("log")
# plt.xscale("log")
# plt.colorbar(hist1, orientation='vertical' )
# plt.colorbar(hist2, orientation='vertical')

# plt.subplot(3,3,5)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.colorbar(hist1, orientation='vertical')
# plt.colorbar(hist2, orientation='vertical')
# plt.yscale("log")
# plt.xscale("log")

# plt.subplot(3,3,6)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.yscale("log")
# plt.xscale("log")
# plt.colorbar(hist1, orientation='vertical', label="all" )
# plt.colorbar(hist2, orientation='vertical', label="undergoing")

# plt.subplot(3,3,7)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.colorbar(hist1, orientation='vertical', )
# plt.colorbar(hist2, orientation='vertical')
# plt.yscale("log")
# plt.xscale("log")

# plt.subplot(3,3,8)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.colorbar(hist1, orientation='vertical', )
# plt.colorbar(hist2, orientation='vertical')
# plt.yscale("log")
# plt.xscale("log")

# plt.subplot(3,3,9)
# hist1 = plt.hexbin(x1, y1, gridsize=30, alpha=1, cmap='Reds', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# hist2 = plt.hexbin(x2, y2, gridsize=30, alpha=0.5, cmap='Blues', xscale="log", yscale="log", bins="log", extent=(0,2, 0,2))
# plt.colorbar(hist1, orientation='vertical', label="all" )
# plt.colorbar(hist2, orientation='vertical', label="undergoing")
# plt.yscale("log")
# plt.xscale("log")
# # plt.margins(0.1) # Uncomment this if hex bins are partially outside of plot limits

# plt.subplots_adjust(top=0.93, bottom=0.07, left=0.1, right=0.9, hspace=0.3, wspace=0.3)

# plt.show()
xxx = np.arange(0.1, 1, 0.1)
xx = np.arange(0.1, 1e2, 1)
print(xx)
print()
yy = np.arange(0.1, 1e2, 1)
hist2 = plt.hexbin(np.concatenate([xxx,xx]), np.concatenate([yy, xxx]), gridsize=30, cmap='Blues', xscale="log", yscale="log", bins="log")
plt.colorbar(hist2, orientation='vertical', label="all" )
# plt.yscale("log")
# plt.xscale("log")
plt.show()

# Generate some data
x = np.random.randn(1000)
y = np.random.randn(1000)

# Create the hexbin plot
fig, ax = plt.subplots()
ax.hexbin(x, y, gridsize=20, mincnt=1)

# Add marginal distributions
ax.hist(x, bins=20, density=True, color='blue', alpha=0.5, )
ax.hist(y, bins=20, density=True, orientation='horizontal', color='red', alpha=0.5)

# Add a title and labels
ax.set_title('Hexbin Plot with Marginal Distributions')
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Show the plot
plt.show()