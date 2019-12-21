The files within ./Identify_Horizon_Python and ./Estimate_Pose_MATLAB
contains all of the data and code necessary to run the pose estimation
algorithm. 

In particular, Algorithm 1 on page 7 of ./../DiffGeo_Final.pdf is
implemented in lines 152-260 of ./Estimate_Pose_MATLAB/horiz_filter.m. This is a MATLAB file
and its required data has been stored in ./Estimate_Pose_MATLAB/data_for_horiz_filter.mat.

The horizon estimation algorithm in Appendix A.2, page 11 of
./../DiffGeo_Final.pdf is implemented in
./Identify_Horizon_Python/process_rocket_vids.py.