# Script to rsync summary files from DAGAGE2 to BP1
# Run this before process_gcwerks, to get the most up-to-date data

# AGAGE GC
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/data/* /work/chxmr/shared/obs_raw/AGAGE_GCWerks/data/
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/data-gcms/* /work/chxmr/shared/obs_raw/AGAGE_GCWerks/data-gcms/
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/air-history/baseline_flags/* /work/chxmr/shared/obs_raw/AGAGE_GCWerks/baseline_flags/

# DECC/GAUGE/AGAGE CRDS
folders=`ssh dagage2.chm.bris.ac.uk "ls /agage/"`
sites="barbados bilsdale heathfield ridgehill tacolneston"

for folder in $folders
do
    for site in $sites
    do
        if [[ $folder == *"$site-picarro"* ]]; then
            rsync -avh dagage2.chm.bris.ac.uk:/agage/${folder}/results-gcwerks/*.nc /work/chxmr/shared/obs_raw/AGAGE_GCWerks/data-picarro/$site/
        fi
    done
done

# Copy ICOS data
rsync -avh dagage2.chm.bris.ac.uk:/agage/macehead-picarro/results-icos/* /work/chxmr/shared/obs_raw/AGAGE_GCWerks/data-icos/macehead/
rsync -avh dagage2.chm.bris.ac.uk:/agage/angus-picarro/results-icos/* /work/chxmr/shared/obs_raw/AGAGE_GCWerks/data-icos/angus/