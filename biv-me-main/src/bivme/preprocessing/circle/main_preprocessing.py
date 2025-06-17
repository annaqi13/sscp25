import os
import numpy as np
import csv
from pathlib import Path

from bivme.preprocessing.circle.utils import *
from bivme.fitting.GPDataSet import GPDataSet


def do_preprocessing(folder, initial_gpfile, initial_sliceinfo, **kwargs):
    """
    Author: Laura Dal Toso
    Date: 20/10/2022
    --------------------------------------------------------------
    This is the main function, that does the pre-processing of GPFiles and SliceInfoFiles to make them
    compatible with the BiVFitting scripts
    -------------------------------------------------------------

    Input:
        - folder where the GPFile and SliceInfoFile are stored

    Output:
        - processed GPFile and SliceInfoFile, ready for the fitting

    """

    if "iter_num" in kwargs:
        iter_num = kwargs.get("iter_num", None)
        pid = os.getpid()
        # assign a new process ID and a new CPU to the child process
        # iter_num corresponds to the id number of the CPU where the process will be run
        os.system("taskset -cp %d %d" % (iter_num, pid))

    if "id_Frame" in kwargs:
        # acquire .csv file containing patient_id, ES frame number, ED frame number if present
        case_frame_dict = kwargs.get("id_Frame", None)

    # First check of SliceInfo and GPFile structure:
    gpfile, sliceinfofile = ReformatFiles(
        folder, initial_gpfile, initial_sliceinfo, temporal_matching
    )

    # chose which frames to upload from the GPFile
    all_frames = pd.read_csv(os.path.join(folder, gpfile), sep="\t")
    frames_to_fit = sorted(np.unique([i[6] for i in all_frames.values]))
    case = os.path.basename(os.path.normpath(folder))

    if do_landmarks_tracking == True:
        print("landmark checking ..")

        data_set = tracking_landmarks(
            folder, gpfile, sliceinfofile, output_csv=landmarks_csv
        )
        if clean_contours == True:
            print("Cleaning contours...", end="")
            final_data_set = Clean_contours(folder, data_set, "GPFile_clean.txt")
            print("done")

    elif do_landmarks_tracking == False and clean_contours == True:
        print("Cleaning contours...", end="")
        filename = os.path.join(folder, gpfile)
        filenameInfo = os.path.join(folder, sliceinfofile)

        data_set = []
        for num in frames_to_fit:
            data_set.append(
                (
                    num,
                    GPDataSet(
                        filename, filenameInfo, case, sampling=1, time_frame_number=num
                    ),
                )
            )

        final_data_set = Clean_contours(folder, data_set, "GPFile_clean.txt")
        print("done")

    else:
        final_data_set = []
        for num in frames_to_fit:
            final_data_set.append(
                (
                    num,
                    GPDataSet(
                        os.path.join(folder, gpfile),
                        os.path.join(folder, sliceinfofile),
                        case,
                        sampling=1,
                        time_frame_number=num,
                    ),
                )
            )

    if find_EDES_frames == True:
        print("Finding ES frame ..")
        findED_ESframe(case, final_data_set)

    DoneFile = Path(os.path.join(folder, "Done.txt"))
    DoneFile.touch(exist_ok=True)

    return DoneFile


if __name__ == "__main__":
    # set directory containing GPFile and SliceInfoFile
    dir_gp = r"R:\resmed201900006-biomechanics-in-heart-disease\Sandboxes\Josh\bivme\test\gpfiles"

    # set list of cases to process
    caselist = ["SCMR_2_corrected", "SCMR_3", "SCMR_4", "SCMR_5"]
    casedirs = [Path(dir_gp, case).as_posix() for case in caselist]

    # check that all landmarks are present at each frame?
    do_landmarks_tracking = False

    # clean up contours (delete basal points)?
    clean_contours = True

    # find which frames are ED and ES?
    find_EDES_frames = False

    # resample contours due to different temporal resolutions between slices?
    temporal_matching = False
    
    # do you want to add valve points from CIM?
    add_cim_valve_points = False
    
    # do you want to add manual slice corrections from CIM?
    add_cim_slice_corrections = False

    initial_gpfile = "GPFile.txt"
    initial_sliceinfo = "SliceInfoFile.txt"

    if do_landmarks_tracking == True:
        fieldnames = [
            "patient",
            "frames",
            "MITRAL_VALVE",
            "TRICUSPID_VALVE",
            "AORTA_VALVE",
            "APEX_POINT",
        ]
        landmarks_csv = "./results/landmarks.csv"

        with open(landmarks_csv, "w") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

    if find_EDES_frames == True:
        with open("./results/case_id_frame.csv", "w") as f:
            writer = csv.DictWriter(f, fieldnames=["patient", "ED", " ES measured"])
            writer.writeheader()

    # start processing...
    [do_preprocessing(folder, initial_gpfile, initial_sliceinfo) for folder in casedirs]

    # if you want to correct GP files using CIM RVLV
    # Add CIM valve points and slice corrections

    if add_cim_valve_points == True or add_cim_slice_corrections == True:
        cim_data = r"R:\resmed201900006-biomechanics-in-heart-disease\Sandboxes\Josh\bivme\test\cim"
        image_ptrs = r"R:\resmed201900006-biomechanics-in-heart-disease\Sandboxes\Josh\bivme\test\images"
        cim_offsets = r"R:\resmed201900006-biomechanics-in-heart-disease\Sandboxes\Josh\biv\chicago\cim_offsets" # TODO: don't need this anymore since we can extract offsets from the cim_data folder
        initial_gpfile = 'GPFile_proc.txt'
        initial_sliceinfo = 'SliceInfoFile_proc.txt'

        [CIM_Correction(folder, initial_gpfile, initial_sliceinfo, cim_data,image_ptrs,cim_offsets,add_cim_valve_points, add_cim_slice_corrections) for folder in casedirs]

        # if you want to clean the contours AFTER adding CIM valve points
        clean_contours = False
        initial_gpfile = 'GPFile_cim.txt'
        initial_sliceinfo = 'SliceInfoFile_cim.txt'
        
        if clean_contours == True:
            [clean_CIM(folder, initial_gpfile, initial_sliceinfo) for folder in casedirs]


    
