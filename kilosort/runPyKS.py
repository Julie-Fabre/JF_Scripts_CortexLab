from pathlib import Path
from pykilosort import run, add_default_handler, np1_probe, np2_probe
from pykilosort.io.probes import Probe
import sys
from mtscomp import decompress
import os
import scipy as sp
import numpy as np

add_default_handler(level='INFO') # print output as the algorithm runs
#probe_dir = {'1':np1_probe(), '2':np2_probe()}


mouse_date = [
                ["JF070","2022-06-11"],
                ["JF070","2022-06-12"],
                ]

RUN = True# False
keepKS2 = True
pykilosortMe = RUN

def custom_probe(xc, yc, chan_tot = 385):
    """ Returns a Neuropixels probe as a Bunch object for use in pykilosort """
    probe_args = {
        'n_channels_total': chan_tot,
        'channel_map': np.arange(384),
        'xcoords': xc,
        'ycoords': yc,
        'channel_groups': np.zeros(384),
        'sample_shifts': np.tile(np.repeat(np.arange(12)/12, 2), 16),
    }

    return Probe(**probe_args)


# define mouse paths
server_path = "/home/netshare/zinu" # change me
mouse_paths = [Path(server_path) / f"{mouse}/{date}/ephys" for mouse, date in mouse_date if len(mouse)>0] # change my formating

# load channel maps
chanmap_path = Path("/home/julie/Dropbox/MATLAB/onPaths/JF_Scripts_CortexLab/kilosort/channelMaps") # change me
chanmap_path_mat = [f.name for f in chanmap_path.iterdir() if ("hanMap" in f.name and ".mat" in f.name)]
channel_maps = {}
for mat in chanmap_path_mat:
    try:
        mat_file = sp.io.loadmat(chanmap_path / mat, mdict=None, appendmat=True)
    except:
        raise ValueError(f"Impossible to load {chanmap_path / mat} for some reason.")

    if len(mat_file['xcoords']) == 385: # remove flipper (sync) channel
        channel_maps[mat] = {
        "xcoords":mat_file['xcoords'].squeeze()[:-1],
        "ycoords":mat_file['ycoords'].squeeze()[:-1],
        }
    else: # flipper not recorded in old open ephys 3a data
        channel_maps[mat] = {
        "xcoords":mat_file['xcoords'].squeeze(),
        "ycoords":mat_file['ycoords'].squeeze(),
        }

# load meta files
def read_metafile(datapath):

    metafile = [f for f in Path(datapath).iterdir() if ".ap.meta" in f.name and "tcat" not in f.name]
    #   isOpenEphys3A = False

    # if len(metafile) == 1: #check if data is in open ephys format
    #     datapath = Path(datapath / "experiment1")
    #     print(f"-datapath: {datapath}.")
    #     if datapath.exists():
    #         isOpenEphys3A = True
    #         print("open ephys 3a")
    #     else:
    #         print("no")

    assert len(metafile) == 1, f".ap.meta file not found or too many .ap.meta at {datapath}!"

    meta_glx = {}

    with open(metafile[0], 'r') as f:
        for ln in f.readlines():
            tmp = ln.split('=')
            k, val = tmp[0], ''.join(tmp[1:])
            k = k.strip()
            val = val.strip('\r\n')
            if '~' in k:
                meta_glx[k] = val.strip('(').strip(')').split(')(')
            else:
                try:  # is it numeric?
                    meta_glx[k] = float(val)
                except:
                    meta_glx[k] = val

        return meta_glx

# iterate over mouse paths
for p in mouse_paths:
    # iterate over site paths
            # find npy files and skip if any found
    if p.exists() == False:
        print(f"No ephys files at {p}.")
        continue

    for subdir in p.iterdir():
        subdir = subdir.name
        pykilosortMe = RUN
        # find whether subdir is a site
        if "-shank" in str(subdir): continue #histology folder
        if "site" not in str(subdir): continue #no ephys, or improperly named
        input_path = p / subdir
        print(f"\n- input path: {input_path}")

        # Find channel maps
        openEphys_path = Path(input_path / "experiment1"/ "recording1" / "continuous" / "Neuropix-3a-100.0")
        print(f"-open ephys path: {openEphys_path}.")
        if openEphys_path.exists():
            isOpenEphys3A = True
            assert('pykilosort can''t handle open ephys 3A data')
        else:
            isOpenEphys3A = False

        if not isOpenEphys3A:
            meta_glx = read_metafile(input_path)
            imro_name = meta_glx["imRoFile"]
        else:
            imro_name = "NPtype3A"

        if not imro_name:
            print('no imro file listed in meta file, default was used')
            prb_type = meta_glx["imDatPrb_type"]
            if prb_type == 0:
                imro_name = "NPtype3B_"
            elif prb_type == 21:
                imro_name = "NPtype21_"
            elif prb_type == 24:
                imro_name = "NPtype24_shank0"
            else:
                assert('no imro file defined')


        if "NPtype24_hStripe_" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype21_" in imro_name:
            chanmap = "chanMapNP2_1Shank_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype24_hStripe_shanks01" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_shank01_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype24_hStripe_shanks23" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_shank23_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype24_shank0" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_shank0_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype24_shank1" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_shank1_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype24_shank2" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_shank2_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype24_shank3" in imro_name:
            chanmap = "chanMapNP2_4Shank_bottRow_shank3_flipper0x2Emat_kilosortChanMap.mat"
        elif "NPtype3B_" in imro_name:
            chanmap = "chanMapNP1_bottRow_flipper.mat"
        elif "NPtype3A" in imro_name:
            chanmap = "chanMapNP1_bottRow.mat"
            input_path = openEphys_path
        else:
            print(f"No pattern found in {imro_name}!")


        xc = channel_maps[chanmap]["xcoords"]
        yc = channel_maps[chanmap]["ycoords"]
        if not isOpenEphys3A:
            chan_tot = meta_glx["nSavedChans"]
        else:
            chan_tot = 384
        probe = custom_probe(xc, yc, chan_tot)

        # find probe version based on binary name
        # prb_version = {False:'1', True:'2'}[".imec" in binary]
        # probe = probe_dir[prb_version]
        # print(f"- probe version: {prb_version}")


        # define output path
        output_path_old = p / "kilosort2" / subdir
        #print(f"- output path: {output_path1}")

        # find npy files and skip if any found
        if output_path_old.exists() and keepKS2:
            np_files = [f for f in output_path_old.iterdir() if ".npy" in str(f)]
            if len(np_files) > 0:
                print(f"Dataset at {output_path_old} already spike sorted with kilosort 2.")
                pykilosortMe = False;

        # define output path
        output_path = p / "pykilosort" / subdir
        output_path_afterKS = p / "pykilosort" / subdir / "output"
        print(f"- output path: {output_path}")

        # find npy files and skip if any found
        if output_path.exists():
            np_files = [f for f in output_path_afterKS.iterdir() if ".npy" in str(f)]
            if len(np_files) > 0:
                print(f"Dataset at {output_path} already spike sorted.")
                pykilosortMe = False;

        # find binary file
        binaries = [f for f in input_path.iterdir() if ((".ap.bin" in str(f) or ".ap.cbin" in str(f) or "continuous.dat" in str(f)) and "tcat" not in str(f))]
        assert len(binaries) != 0, f"No file was found at {input_path}!"
        assert len(binaries) == 1, f"More than one binary file found at {input_path}!"
        binary = str(binaries[0])
        print(f"- binary file: {binary}")

        # decompress data, save on drive (extract sync with matlab scripts later)
        if '.cbin' in binary:
            sync_path = output_path  / "sync.mat"
            sync_path_old = output_path_old  / "sync.mat"
            print(f"- sync path: {sync_path}")
            if sync_path.exists() or sync_path_old.exists() :
                print(f"- sync already extracted")
                final_input_path = input_path / binary
            else:
                binary_ch = binary[:-5]
                binary_ch = binary_ch + '.ch'
                local_save_path = os.path.normpath(p)
                local_save_path_split = local_save_path.split(os.sep)
                local_bin_path = os.path.normpath(binary)
                local_bin_path_split = local_bin_path.split(os.sep)
                local_bin_path_name = local_bin_path_split[-1][:-5]
                local_bin_path_name = local_bin_path_name + '.bin'
                local_decompressed_data_path = Path('/media/julie/Elements/' + local_save_path_split[-3] + '/' + local_save_path_split[-2] + '/' + local_save_path_split[-1] + '/' + local_bin_path_name)
                print(f"- decompressed save path: {local_decompressed_data_path}")
                if local_decompressed_data_path.exists():
                    print(f"Data already decompressed at: {local_decompressed_data_path}")
                    final_input_path = local_decompressed_data_path
                else:
                    if RUN:
                        saveDecompPath = Path('/media/julie/Elements/' + local_save_path_split[-3] + '/' + local_save_path_split[-2] + '/' + local_save_path_split[-1]);
                        if saveDecompPath.exists():
                            print('directory aleady there')
                        else:
                            os.makedirs(Path('/media/julie/Elements/' + local_save_path_split[-3] + '/' + local_save_path_split[-2] + '/' + local_save_path_split[-1]))

                        decompress(binary,cmeta=binary_ch,out=local_decompressed_data_path)
                        final_input_path = local_decompressed_data_path
                    else:
                        print(f"Skipping data decompression and kilosorting of {input_path}")
        else:
            final_input_path = input_path / binary

        # run pykilosort
        if RUN and pykilosortMe:
            try:
                print(final_input_path)
                run(final_input_path, dir_path=output_path, probe=probe)
            except:
                print(f"Error kilosorting data at {output_path}.")
