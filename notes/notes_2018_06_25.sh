# notes for looking for upper and lower discontinuities in phased files

# find line number of de novo
grep -n "de novo end chromosome position" filename (from bed)

# looking through lines
sed -n "start,(end)p;(end+1)q" file_name | cut -f10 | grep -E "0/1|1/0"

# once you have result, paste into quotes of
grep -n "" filename
