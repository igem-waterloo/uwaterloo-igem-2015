import errno
import os
import re
import shutil
import subprocess


def natural_sort(unsorted_list):
    """Used to sort lists like a human
    e.g. [1, 10, 11, 2, 3] to [1, 2, 3, 10, 11]
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(unsorted_list, key=alphanum_key)


def copy_and_rename_plots(plot_lattice_dir, output_dir):
    """
    Given a plot_lattice or plot_grid directory, copies it and renames the files in an order that
    ffmpeg.exe likes, e.g. (0001, 0002, 0003, ... , 0010, 0011, ...)
    Notes:
        - assumes less than 10000 files are being copied (for ffmpeg simplicity)
    """
    # copy the folder
    try:
        shutil.copytree(plot_lattice_dir, output_dir)
    except OSError as exc:
        if exc.errno == errno.ENOTDIR:
            shutil.copy(plot_lattice_dir, output_dir)
        else:
            raise
    # naturally sort the copied files
    unsorted_files = os.listdir(plot_lattice_dir)
    assert len(unsorted_files) <= 99999  # assume <= 5 digits later
    sorted_files = natural_sort(unsorted_files)
    # rename the files accordingly
    basename = "genome_"
    filetype = ".png"
    for i, filename in enumerate(sorted_files):
        num = "%05d" % i
        newname = basename + num + filetype
        os.rename(os.path.join(output_dir, filename), os.path.join(output_dir, newname))
    return


def make_video_ffmpeg(plot_genome_dir, output_path, fps=15, rename_flag=False, ffmpeg_dir=None):
    """Makes a video using ffmpeg - also copies the lattice plot dir, changes filenames, and deletes the copy
    Args:
        plot_genome_dir: source directory
        output_path: path and filename of the output video
        fps: frames per second of the output video
        rename_flag: [default: False] rename to format that ffmpeg likes
        ffmpeg_dir: [default: None] location of the ffmpeg directory (root where bin containing ffmpeg.exe is)
    Returns:
        None
    Notes:
        - assumes ffmpeg has been extracted on your system and added to the path
        - if it's not added to path, point to it (the directory containing ffmpeg bin) using ffmpeg_dir arg
        - assumes less than 10000 images are being joined (for ffmpeg simplicity)
    """
    source_dir = plot_genome_dir
    if rename_flag:  # make temp dir
        temp_plot_dir = os.path.join(plot_genome_dir, os.pardir, "temp")
        copy_and_rename_plots(plot_genome_dir, temp_plot_dir)
        source_dir = temp_plot_dir

    # make video
    command_line = ["ffmpeg", "-framerate", "%d" % fps, "-i", os.path.join(source_dir, "genome_%05d.png"),
                    "-c:v", "libx264", "-r", "%d" % fps, "-pix_fmt", "yuv420p", "%s" % output_path]
    if ffmpeg_dir is not None:
        app_path = os.path.join(ffmpeg_dir, "bin", "ffmpeg.exe")
        sp = subprocess.Popen(command_line, executable=app_path, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    else:
        sp = subprocess.Popen(command_line, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = sp.communicate()
    print out, err, sp.returncode

    # delete temp dir
    if rename_flag:
        shutil.rmtree(temp_plot_dir)

    return
