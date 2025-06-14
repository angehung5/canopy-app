"""
Created on Sun Jun 4  2023
Updated on Mon Oct 16 2023: Use daily gfs.canopy files
Updated on Tue Nov 7  2023: Enable multiple times as user argument
Updated on Tue Apr 2  2024: Remove wget functions, all data must be from local files
Updated on Fri May 31 2024: Replace scipy griddata with monet (pyresample)
Updated on Mon Feb 24 2025: Enable GFSv16 downloading from AWS; Add "ozone_w126"
Updated on Tue Jun 10 2025: Enable to specify sources of each variable

Author: Wei-Ting Hung
"""

import os
import subprocess
import sys
from datetime import datetime, timedelta, timezone

import monet  # noqa: F401
import numpy as np
import xarray as xr
from netCDF4 import Dataset
from pysolar.solar import get_altitude

"""User Arguments"""
# Time: yyyymmddhhfff
timelist = sys.argv[1]
timelist = np.array(timelist.split(",")).astype(str)

"""User Options"""
path = "./input"  # work directory
ref_lev = 10  # reference height (m, a.g.l.)
frp_src = 0  # frp data source (0: local fire product; 1: 12 month climatology; 2: all ones when ifcanwaf=.FALSE.)
can_src = 0  # canopy data source (0: pre-generated daily global file based on year 2022; 1: global canopy data file from AWS; 2: user specified)


# ------------------------------ ATTENTION -------------------------------- #
# UPDATE - June 10 2025                                                     #
# If user specified canopy data is used (can_src=2), please specify file    #
# locations in function "find_user_canopy". Data reading and processing     #
# can be specifiedd in function "read_user_canopy". Please check function   #
# "read_aws_canopy" as an example.                                          #
#                                                                           #
# ------------------------------------------------------------------------- #
# If local FRP is used (frp_src=0,1), archived GBBEPx files since 2020 are  #
# available for GMU HOPPER users. For users outside GMU, parameter "f_frp"  #
# and function "read_frp_local" need to be modified accordingly.            #
# Function read_frp_local is designed for reading GBBEPx by default.        #
#                                                                           #
# Recent GBBEPx v3 files (~ < 3 months) are available for download at:      #
# https://www.ospo.noaa.gov/Products/land/gbbepx/                           #
# GBBEPx_all01GRID.emissions_v003_"+YY+MM+DD+".nc"                          #
#                                                                           #
# Archived GBBEPx files since 2020 are available for GMU HOPPER users.      #
#                                                                           #
# 12 month climatological FRP will not be respective of actual conditions.  #
# Only use it when users do not need the limited WAF application for fires. #
# ------------------------------------------------------------------------- #


"""Global Settings"""
# domain
lat_lim = [-90, 90]
lon_lim = [0, 360]


# required variables
loclist = ["grid_xt", "lon", "grid_yt", "lat", "time"]
metlist = [
    "ugrd10m",
    "vgrd10m",
    "fricv",
    "sfcr",
    "vtype",
    "sotyp",
    "pressfc",
    "dswrf",
    "shtfl",
    "tmpsfc",
    "tmp2m",
    "spfh2m",
    "hpbl",
    "prate_ave",
]
canlist = [
    "lai",
    "clu",
    "canfrac",
    "ch",
    "pavd",
    "ozone_w126",
    "mol",
    "csz",
    "frp",
    "href",
]


# constants
fill_value = 9.99e20  # fill value
den = 1.18  # air density (kg/m3)
Cp = 1004  # specific heat capacity of Air at 25C (J/kg/K)
K = 0.4  # Von Karman constant
g = 9.8  # gravitational acceleration (m/s2)


# functions
def read_varatt(var):
    attname = var.ncattrs()
    att = [var.getncattr(X) for X in attname]
    return attname, att


def write_varatt(var, attname, att):
    for X in np.arange(len(attname)):
        if attname[X] == "_FillValue":
            continue
        elif attname[X] == "missing_value":
            value = np.float32(att[X])
            value = np.round(value / 1e15) * 1e15
            var.setncattr(attname[X], value)
        else:
            var.setncattr(attname[X], att[X])


def find_user_canopy(year):
    # Please specify the location of local user canopy data here.
    flist = {
        "lai": "/groups/ESS/whung/Alldata/Global_canopy/grid1km/canopy_leaf_area_index."
        + year
        + ".0.01.nc",
        "clu": "/groups/ESS/whung/Alldata/Global_canopy/grid1km/canopy_clumping_index.2001_2017.0.01.nc",
        "canfrac": "/groups/ESS/whung/Alldata/Global_canopy/grid1km/canopy_green_vegetation_fraction."
        + year
        + ".0.01.nc",
        "ch": "/groups/ESS/whung/Alldata/Global_canopy/grid1km/canopy_height.2020.0.01.nc",
        "pavd": "/groups/ESS/whung/Alldata/Global_canopy/grid1km/canopy_plant_area_volume_density.2019_2023.0.01.nc",
    }
    return flist


def read_gfs_climatology(filename, basefile, varname):
    readin = xr.open_dataset(filename)
    readin = readin.set_coords(["lat", "lon"]).rename(
        {"grid_xt": "x", "grid_yt": "y", "lat": "latitude", "lon": "longitude"}
    )

    nlev = len(readin.lev.data)

    if varname == "pavd":
        DATA = np.empty([nlev, basefile.zc.data.shape[1], basefile.zc.data.shape[2]])
        for ll in np.arange(nlev):
            DATA[ll, :, :] = (
                basefile["zc"].monet.remap_nearest(readin[varname][0, ll, :, :]).data
            )
    else:
        DATA = basefile["zc"].monet.remap_nearest(readin[varname][0, :, :]).data

    readin.close()
    DATA[np.isnan(DATA)] = 0
    DATA[DATA < 0] = 0
    return DATA


def read_aws_canopy(filename, basefile, varname, month):
    readin = xr.open_dataset(filename)
    readin = readin.rename({"jdim": "y", "idim": "x", "jdim_p1": "y_p1"})

    grid_xt, grid_yt = np.meshgrid(readin["lon"].data, readin["lat"].data)
    yt = xr.DataArray(grid_yt, dims=["y", "x"], name="latitude")
    xt = xr.DataArray(grid_xt, dims=["y", "x"], name="longitude")
    readin["latitude"] = yt
    readin["longitude"] = xt
    readin = readin.set_coords(["latitude", "longitude"])

    if varname == "lai":
        DATA = (
            basefile["zc"]
            .monet.remap_nearest(readin["canopy_leaf_area_index"][int(month) - 1, :, :])
            .data
        )
        DATA[DATA < 0] = 0
    elif varname == "clu":
        DATA = (
            basefile["zc"]
            .monet.remap_nearest(readin["canopy_clumping_index"][int(month) - 1, :, :])
            .data
        )
        DATA[DATA < 0] = 0
    elif varname == "canfrac":
        DATA = (
            basefile["zc"]
            .monet.remap_nearest(
                readin["canopy_green_vegetation_fraction"][int(month) - 1, :, :]
            )
            .data
        )
        DATA[DATA < 0] = 0
        DATA[DATA > 1] = 1
    elif varname == "ch":
        DATA = basefile["zc"].monet.remap_nearest(readin["canopy_height"][0, :, :]).data
        DATA[DATA < 0] = 0
    elif varname == "pavd":
        nlev = len(readin.layer_top.data)
        DATA = np.empty([nlev, basefile.zc.data.shape[1], basefile.zc.data.shape[2]])
        for ll in np.arange(nlev):
            DATA[ll, :, :] = (
                basefile["zc"]
                .monet.remap_nearest(
                    readin["canopy_plant_area_volume_density"][0, ll, :, :]
                )
                .data
            )
        DATA[DATA < 0] = 0
    readin.close()
    return DATA


def read_user_canopy():
    # Please specify the necessary processing of user canopy data here.
    # Recommanded data processes include but not limit to: data read in, unit conversion and gridding.
    DATA = np.empty(1)
    return DATA


def read_frp_local(filename, basefile):
    readin = xr.open_dataset(filename)
    readin = readin.rename({"Latitude": "y", "Longitude": "x"})
    readin["x"] = readin["x"].where(readin["x"] > 0, readin["x"] + 360)

    grid_xt, grid_yt = np.meshgrid(readin["x"].data, readin["y"].data)
    yt = xr.DataArray(grid_yt, dims=["y", "x"], name="latitude")
    xt = xr.DataArray(grid_xt, dims=["y", "x"], name="longitude")
    readin["latitude"] = yt
    readin["longitude"] = xt
    readin = readin.set_coords(["latitude", "longitude"])

    DATA = basefile["zc"].monet.remap_nearest(readin["MeanFRP"][0, :, :]).data
    readin.close()
    return DATA


starttime = datetime.now()
print("------------------------------------")
print("---- Global input pre-process start!", starttime.strftime("%Y/%m/%d %H:%M:%S"))
print("------------------------------------")


for inputtime in timelist:
    if (
        (len(inputtime) != 13)
        or (int(inputtime[4:6]) > 12)
        or (int(inputtime[6:8]) > 31)
        or (int(inputtime[8:10]) > 23)
        or (int(inputtime[10:]) > 23)
    ):
        print("---- Unidentified time format: " + inputtime + ". Terminated!")
        exit()
    else:
        print("---- " + inputtime[:10] + " f" + inputtime[10:] + " processing...")
        print("------------------------------------")

    """Settings"""
    # date/time
    YY = inputtime[:4]  # year
    MM = inputtime[4:6]  # month
    DD = inputtime[6:8]  # date
    HH = inputtime[8:10]  # hour in UTC
    FH = inputtime[10:]  # forecast hour, for met file only

    # input/output files
    f_met = (
        path + "/gfs.t" + HH + "z." + YY + MM + DD + ".sfcf" + FH + ".nc"
    )  # gfs met file
    f_can = (
        path
        + "/gfs.canopy.t"
        + HH
        + "z."
        + "2022"
        + MM
        + DD
        + ".sfcf"
        + FH
        + ".global.nc"
    )  # canopy file
    f_output = (
        path + "/gfs.t" + HH + "z." + YY + MM + DD + ".sfcf" + FH + ".canopy.nc"
    )  # output file

    if frp_src == 0:  # local frp file
        if int(YY + MM + DD) <= 20230510:  # version 3
            f_frp = (
                "/groups/ESS/yli74/data/GBBEPx/ORI/GBBEPx_all01GRID.emissions_v003_"
                + YY
                + MM
                + DD
                + ".nc"
            )
        else:  # version 4
            f_frp = (
                "/groups/ESS/yli74/data/GBBEPx/ORI/GBBEPx-all01GRID_v4r0_blend_"
                + YY
                + MM
                + DD
                + ".nc "
            )
    elif frp_src == 1:  # climatological frp
        f_frp = (
            path + "/gfs.canopy.t" + HH + "z." + "2022" + MM + DD + ".sfcf000.global.nc"
        )

    if can_src == 1:  # global canopy data file from AWS
        f_can_list = {
            "lai": path + "/canopy_leaf_area_index." + YY + ".0.01.nc",
            "clu": path + "/canopy_clumping_index.2001_2017.0.01.nc",
            "canfrac": path + "/canopy_green_vegetation_fraction." + YY + ".0.01.nc",
            "ch": path + "/canopy_height.2020.0.01.nc",
            "pavd": path + "/canopy_plant_area_volume_density.2019_2023.0.01.nc",
        }
    elif can_src == 2:  # user specified local canopy data
        f_can_list = find_user_canopy(YY)

    """Data Check"""
    """Program terminates if required files do not exist."""
    print("---- Checking required files...")
    print("------------------------------------")

    # met file
    if os.path.isfile(f_met) is True:
        print("---- Met file found!")
    else:
        subprocess.run(
            [
                "wget",
                "--no-check-certificate",
                "--no-proxy",
                "-O",
                path + "/gfs.t12z." + YY + MM + DD + ".sfcf" + FH + ".nc",
                "https://noaa-oar-arl-nacc-pds.s3.amazonaws.com/inputs/"
                + YY
                + MM
                + DD
                + "/gfs.t12z.sfcf"
                + FH
                + ".nc",
            ]
        )
        if os.path.isfile(f_met) is True:
            os.chmod(f_met, 0o0755)
            print("---- Met file downloaded!")
        else:
            print("---- No available met data. Terminated!")
            exit()

    # can file
    if os.path.isfile(f_can) is True:
        print("---- Canopy file found!")
    else:
        subprocess.run(
            [
                "wget",
                "--no-check-certificate",
                "--no-proxy",
                "-O",
                path
                + "/gfs.canopy.t12z."
                + "2022"
                + MM
                + DD
                + ".sfcf"
                + FH
                + ".global.nc",
                "https://noaa-oar-arl-nacc-pds.s3.amazonaws.com/inputs/geo-files/"
                + "gfs.canopy.t12z."
                + "2022"
                + MM
                + DD
                + ".sfcf000.global.nc",
            ]
        )
        if os.path.isfile(f_can) is True:
            os.chmod(f_can, 0o0755)
            print("---- Canopy file downloaded!")
        else:
            print("---- No available canopy data. Terminated!")
            exit()
    if can_src == 1:  # global files from AWS
        for key in f_can_list.keys():
            f = f_can_list[key]
            print(f)
            subprocess.run(
                [
                    "wget",
                    "--no-check-certificate",
                    "--no-proxy",
                    "-O",
                    f,
                    "https://noaa-oar-arl-nacc-pds.s3.amazonaws.com/inputs/geo-files/"
                    + f[f.rindex("/") + 1 :],
                ]
            )
            if os.path.isfile(f) is True:
                os.chmod(f, 0o0755)
                print("----", f[f.rindex("/") + 1 :], "downloaded!")
            else:
                print("----", f[f.rindex("/") + 1 :], "not available. Terminated!")
                exit()
    elif can_src == 2:  # user specified
        checklist = [
            os.path.isfile(f_can_list["lai"]),
            os.path.isfile(f_can_list["clu"]),
            os.path.isfile(f_can_list["canfrac"]),
            os.path.isfile(f_can_list["ch"]),
            os.path.isfile(f_can_list["pavd"]),
        ]
        if False in checklist:
            print("---- Invalid user specified canopy data found:")
            print("lai   clu   canfrac   ch   pavd")
            print(checklist)
            print("Terminated!")
            exit()

    # frp file
    if frp_src == 0:  # local source
        if os.path.isfile(f_frp) is True:
            os.system("cp " + f_frp + " " + path)
            if int(YY + MM + DD) <= 20230510:
                f_frp = path + "/GBBEPx_all01GRID.emissions_v003_" + YY + MM + DD + ".nc"
            else:
                f_frp = path + "/GBBEPx-all01GRID_v4r0_blend_" + YY + MM + DD + ".nc"
            os.chmod(f_frp, 0o0755)
            print("---- FRP file found!")
        else:
            print("---- No available FRP file. Terminated!")
            exit()

    if frp_src == 1:  # 12 month climatology frp
        if os.path.isfile(f_frp) is True:
            print("---- FRP file found!")
        else:
            print("---- No available FRP data. Terminated!")
            exit()

        print("-----------!!!WARNING!!!-------------")
        print("---!!!Climatological FRP is used!!!--")

    os.system("cp " + f_met + " " + f_output)  # copy gfs met file for appending

    """Reading dimensions"""
    print("------------------------------------")
    print("---- Checking variable dimensions...")
    print("------------------------------------")
    basefile = xr.open_dataset(f_met)
    basefile = basefile.set_coords(["lat", "lon"]).rename(
        {"grid_xt": "x", "grid_yt": "y", "lat": "latitude", "lon": "longitude"}
    )

    # dimension sizes
    ntime = len(basefile["time"].data)
    nlat = len(basefile["y"].data)
    nlon = len(basefile["x"].data)

    # var check
    print("time", basefile["time"].data.shape)
    print("grid_yt", basefile["y"].data.shape)
    print("grid_xt", basefile["x"].data.shape)
    print("lat", basefile["latitude"].data.shape)
    print("lon", basefile["longitude"].data.shape)

    """Adding canvar"""
    print("------------------------------------")
    print("---- Generating canopy variables...")
    print("------------------------------------")

    for i in np.arange(len(canlist)):
        varname = canlist[i]

        print("---- " + varname + " processing...")

        if varname == "lai":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Leaf area index", "m^2/m^2", fill_value]
            if can_src == 0:
                DATA = read_gfs_climatology(f_can, basefile, "lai")
            elif can_src == 1:
                DATA = read_aws_canopy(f_can_list["lai"], basefile, "lai", MM)
            elif can_src == 2:
                DATA = read_user_canopy()

        elif varname == "clu":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Canopy clumping index", "none", fill_value]
            if can_src == 0:
                DATA = read_gfs_climatology(f_can, basefile, "clu")
            elif can_src == 1:
                DATA = read_aws_canopy(f_can_list["clu"], basefile, "clu", MM)
            elif can_src == 2:
                DATA = read_user_canopy()

        elif varname == "canfrac":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Forest fraction of grid cell", "none", fill_value]
            if can_src == 0:
                DATA = read_gfs_climatology(f_can, basefile, "canfrac")
            elif can_src == 1:
                DATA = read_aws_canopy(f_can_list["canfrac"], basefile, "canfrac", MM)
            elif can_src == 2:
                DATA = read_user_canopy()

        elif varname == "ch":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Canopy height above the surface", "m", fill_value]
            if can_src == 0:
                DATA = read_gfs_climatology(f_can, basefile, "ch")
            elif can_src == 1:
                DATA = read_aws_canopy(f_can_list["ch"], basefile, "ch", MM)
            elif can_src == 2:
                DATA = read_user_canopy()

        elif varname == "pavd":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Plant area volume density profile", "m2/m3", fill_value]
            if can_src == 0:
                DATA = read_gfs_climatology(f_can, basefile, "pavd")
            elif can_src == 1:
                DATA = read_aws_canopy(f_can_list["pavd"], basefile, "pavd", MM)
            elif can_src == 2:
                DATA = read_user_canopy()

        elif varname == "ozone_w126":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Ozone W126 index", "ppm-hours", fill_value]
            DATA = read_gfs_climatology(f_can, basefile, "ozone_w126")

        elif varname == "mol":
            # Reference:
            # Essa 1999, ESTIMATION OF MONIN-OBUKHOV LENGTH USING RICHARDSON AND BULK RICHARDSON NUMBER
            # https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/118/37118528.pdf
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Monin-Obukhov length", "m", fill_value]

            readin = Dataset(f_met)
            t2m = np.squeeze(readin["tmp2m"][:])
            fricv = np.squeeze(readin["fricv"][:])
            shtfl = np.squeeze(readin["shtfl"][:])

            DATA = (-1 * den * Cp * t2m * (fricv**3)) / (K * g * shtfl)
            DATA[DATA > 500] = 500
            DATA[DATA < -500] = -500

            del [readin, t2m, fricv, shtfl]

        elif varname == "csz":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Cosine of solar zenith angle", "none", fill_value]

            lat = basefile.latitude.data
            lon = basefile.longitude.data

            time_conv = datetime(
                int(YY), int(MM), int(DD), int(HH), 0, 0, 0, tzinfo=timezone.utc
            ) + timedelta(hours=int(FH))
            sza = 90 - get_altitude(lat, lon, time_conv)
            DATA = np.cos(sza * 0.0174532925)  # degree to radian

            del [lat, lon, time_conv, sza]

        elif varname == "frp":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Mean fire radiative power", "MW", fill_value]

            if frp_src == 1:  # 12 month climatology
                DATA = read_gfs_climatology(f_can, basefile, "frp")
            elif frp_src == 2:  # ifcanwaf=.FALSE.
                DATA = np.empty(lat.shape)
                DATA[:] = 1
            else:
                DATA = read_frp_local(f_frp, basefile)

        elif varname == "href":
            ATTNAME = ["long_name", "units", "missing_value"]
            ATT = ["Reference height above the surface", "m", fill_value]
            DATA = np.empty([nlat, nlon])
            DATA[:] = ref_lev

        # var check
        print("Dimension/Attributes:")
        print(DATA.shape)
        print(ATTNAME)
        print(ATT)

        # adding to output file
        if varname == "pavd":
            output = Dataset(f_output, "a")
            output.createDimension("level", DATA.shape[0])

            var_lev = output.createVariable("lev", "float", ("level",))
            var_bot = output.createVariable("layer_bottom", "i4", ("level",))
            var_top = output.createVariable("layer_top", "i4", ("level",))
            var = output.createVariable(
                varname,
                "float",
                ("time", "level", "grid_yt", "grid_xt"),
                fill_value=fill_value,
            )

            write_varatt(var, ATTNAME, ATT)
            write_varatt(
                var_lev,
                ["long_name", "units"],
                ["height of the layer midpoint above the ground", "m"],
            )
            write_varatt(
                var_bot,
                ["long_name", "units"],
                ["height of the layer bottom above the ground", "m"],
            )
            write_varatt(
                var_top,
                ["long_name", "units"],
                ["height of the layer top above the ground", "m"],
            )

            var_lev[:] = np.arange(2.5, 70, 5)
            var_bot[:] = np.arange(0, 65 + 1, 5)
            var_top[:] = np.arange(5, 70 + 1, 5)
            var[:] = DATA

            output.close()
            del [var_bot, var_top]

        else:
            output = Dataset(f_output, "a")

            var = output.createVariable(
                varname, "float", ("time", "grid_yt", "grid_xt"), fill_value=fill_value
            )
            write_varatt(var, ATTNAME, ATT)
            var[:] = DATA

            output.close()

        print("---- " + varname + " complete!")

        del [output, var]
        del [varname, DATA, ATTNAME, ATT]

    print("------------------------------------")
    print("---- " + inputtime + " complete!")
    print("------------------------------------")


endtime = datetime.now()


print("---- Global input pre-process complete!", endtime.strftime("%Y/%m/%d %H:%M:%S"))
print("---- Process time:", str(endtime - starttime))
print("------------------------------------")
