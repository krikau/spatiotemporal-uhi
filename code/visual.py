import matplotlib.pyplot as plt
import cartopy


def temperature_gradients_on_map(temperature_gradients, lon_range, lat_range):
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=cartopy.crs.PlateCarree())
    im = ax.imshow(
        temperature_gradients,
        extent=lon_range + lat_range,
        transform=cartopy.crs.PlateCarree(),
        cmap="coolwarm",
        vmin=-10,
        vmax=10,
    )
    ax.coastlines(resolution="10m")
    ax.set_title("Temperature Gradients on Map")
    plt.colorbar(im, label="Temperature Gradient (°C)")
    plt.show()


def temperature_gradients(temperature_gradients):
    plt.imshow(temperature_gradients, cmap="coolwarm", vmin=-10, vmax=10)
    plt.colorbar(label="Temperature Gradient (°C)")
    plt.title("Temperature Gradients")
    plt.xlabel("Column Index")
    plt.ylabel("Row Index")
    plt.show()
