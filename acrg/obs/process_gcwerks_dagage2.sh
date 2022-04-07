# Script to rsync summary files from DAGAGE2 to BP1
# Run this before process_gcwerks, to get the most up-to-date data

# AGAGE GC
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/data/* /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data/
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/data-gcms/* /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data-gcms/
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/air-history/baseline_flags/* /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/baseline_flags/

# DECC GCMD
rsync -avh dagage2.chm.bris.ac.uk:/agage/summary/gccompare-decc/data/md/* /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data/

# DECC/GAUGE/AGAGE CRDS
folders=`ssh dagage2.chm.bris.ac.uk "ls /agage/"`
sites="barbados bilsdale heathfield ridgehill tacolneston"

for folder in $folders
do
    for site in $sites
    do
        if [[ $folder == *"$site-picarro"* ]]; then
            rsync -avh dagage2.chm.bris.ac.uk:/agage/${folder}/results-gcwerks/*.dat /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data-picarro/$site/
        fi
    done
done

# Get HFD 5310
rsync -avh dagage2.chm.bris.ac.uk:/agage/heathfield-5310/results-gcwerks/*.dat /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data-picarro/heathfield/

# Get TAC LGR
rsync -avh dagage2.chm.bris.ac.uk:/agage/tacolneston-lgr/results-gcwerks/*.dat /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data-picarro/tacolneston/

# Copy ICOS data
rsync -avh dagage2.chm.bris.ac.uk:/agage/macehead-picarro/results-icos/* /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data-icos/macehead/
rsync -avh dagage2.chm.bris.ac.uk:/agage/angus-picarro/results-icos/* /group/chemistry/acrg/obs_raw/AGAGE_GCWerks/data-icos/angus/