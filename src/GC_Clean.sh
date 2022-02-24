# cd ~/Documents/GitHub/GrowthCurveAnalysis/GC_Data

# Need gnumeric package `brew install gnumeric` to do the automated conversion from xls to csv
for j in `ls OGExportDataFiles/xlsx/* | awk -F '[\/.]' '{print $3}'`
	do
	ssconvert OGExportDataFiles/xlsx/$j.xl* OGExportDataFiles/csv/$j.csv
	done

# Lines I used to take the original export files and extract only the data tables with Time, Well, and OD600 data. Keeping them here for posterity. 

sed -n "59,108p"  OGExportDataFiles/csv/GrowthCurve24hr_Ctr9Cdc73_201012_001.csv | awk -F '[,]' '{ $1=""; $3=""; OFS = "," ; print $0 }' | sed 's/,,/,/g' | cut -c 2- > DataTables/CtCdGCData_001.csv

sed -n "59,108p"  OGExportDataFiles/csv/GrowthCurve24h_Paf1Rtf1_20201020_002.csv | awk -F '[,]' '{ $1=""; $3=""; OFS = "," ; print $0 }' | sed 's/,,/,/g' | cut -c 2- > DataTables/PaRtGCData_002.csv

sed -n "59,108p"  OGExportDataFiles/csv/GrowthCurve24h_Paf1Rtf1_20201026_003.csv | awk -F '[,]' '{ $1=""; $3=""; OFS = "," ; print $0 }' | sed 's/,,/,/g' | cut -c 2- > DataTables/PaRtGCData_003.csv

sed -n "59,108p"  OGExportDataFiles/csv/GrowthCurve24hr_20201125_CtLe.csv | awk -F '[,]' '{ $1=""; $3=""; OFS = "," ; print $0 }' | sed 's/,,/,/g' | cut -c 2- > DataTables/CtLeGCData_004.csv

sed -n "59,108p"  OGExportDataFiles/csv/GrowthCurve24h_RtRth_20201203.csv | awk -F '[,]' '{ $1=""; $3=""; OFS = "," ; print $0 }' | sed 's/,,/,/g' | cut -c 2- > DataTables/RtRthGCData_005.csv


sed -n "59,108p"  OGExportDataFiles/csv/GrowthCurve24h_RthLe_20201208.csv | awk -F '[,]' '{ $1=""; $3=""; OFS = "," ; print $0 }' | sed 's/,,/,/g' | cut -c 2- > DataTables/RthLeGCData_006.csv
