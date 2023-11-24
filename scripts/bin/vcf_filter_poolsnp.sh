#!bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH --cpus-per-task 1
#SBATCH --time 10:00:00

input_file=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic/EntNic.vcf.gz
output_file=/faststorage/project/EcoGenetics/BACKUP/population_genetics/collembola/Entomobrya_nicoleti/EntNic/EntNic_filtered.vcf.gz

if [ "${input_file%.*}" == ".gz" ]; then
	input_file="<(zcat "$input_file")"
fi

awk \
'BEGIN{FS=OFS="\t"}
{
len=length($5) 
if ($5 !=  "." && len == 1)
	{for(i=1;i<=NF;i++)
		{if (substr($i,1,3)=="0/1" && substr($i,length($i)-2,3)=="0.0") 
			{sub("0/1", "0/0", $i)}
		}; print
	}
else
	{if ($5 !=  "." && len == 3)
		{for(i=1;i<=NF;i++)
			{if (substr($i,1,4)=="0/1:" && substr($i,length($i)-2,3)=="0.0")
				{sub("0/1", "0/0", $i)}
			if (substr($i,1,4)=="0/2:" && substr($i,length($i)-2,3)=="0.0")
				{sub("0/2", "0/0", $i)}
			if (substr($i,1,6)=="0/1/2:" && substr($i,length($i)-6,7)=="0.0,0.0")
				{sub("0/1/2", "0/0", $i)}
			if (substr($i,1,6)=="0/1/2:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
				{sub("0/1/2", "0/1", $i)}
			if (substr($i,1,6)=="0/1/2:" && substr($i,length($i)-7,3)=="0.0" && substr($i,length($i)-2,3)!="0.0")
				{sub("0/1/2", "0/2", $i)}
			}; print
		}
	else
		{if ($5 !=  "." && len == 5)
			{for(i=1;i<=NF;i++)
				{if (substr($i,1,4)=="0/1:"  && substr($i,length($i)-2,3)=="0.0")
					{sub("0/1", "0/0", $i)}
				if (substr($i,1,4)=="0/2:" && substr($i,length($i)-2,3)=="0.0")
					{sub("0/2", "0/0", $i)}
				if (substr($i,1,4)=="0/3:" && substr($i,length($i)-2,3)=="0.0")
					{sub("0/3", "0/0", $i)}
				if (substr($i,1,6)=="0/1/2:" && substr($i,length($i)-6,7)=="0.0,0.0")
					{sub("0/1/2", "0/0", $i)}
				if (substr($i,1,6)=="0/1/2:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
					{sub("0/1/2", "0/1", $i)}
				if (substr($i,1,6)=="0/1/2:"&& substr($i,length($i)-7,3)=="0.0" && substr($i,length($i)-2,3)!="0.0")
					{sub("0/1/2", "0/2", $i)}
				if (substr($i,1,6)=="0/1/3:" && substr($i,length($i)-6,7)=="0.0,0.0")
					{sub("0/1/3", "0/0", $i)}
				if (substr($i,1,6)=="0/1/3:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
					{sub("0/1/3", "0/1", $i)}
				if (substr($i,1,6)=="0/1/3:" && substr($i,length($i)-7,3)=="0.0" && substr($i,length($i)-2,3)!="0.0")
					{sub("0/1/3", "0/3", $i)}
				if (substr($i,1,6)=="0/2/3:" && substr($i,length($i)-6,7)=="0.0,0.0")
					{sub("0/2/3", "0/0", $i)}
				if (substr($i,1,6)=="0/2/3:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
					{sub("0/2/3", "0/2", $i)}
				if (substr($i,1,6)=="0/2/3:" && substr($i,length($i)-7,3)=="0.0" && substr($i,length($i)-2,3)!="0.0")
					{sub("0/2/3", "0/3", $i)}
				if (substr($i,1,8)=="0/1/2/3:" && substr($i,length($i)-10,11)=="0.0,0.0,0.0") 
					{sub("0/1/2/3", "0/0", $i)}
				if (substr($i,1,8)=="0/1/2/3:" && substr($i,length($i)-6,7)=="0.0,0.0" && substr($i,length($i)-10,3)!="0.0") 
					{sub("0/1/2/3", "0/1", $i)}
				if (substr($i,1,8)=="0/1/2/3:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-11,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
					{sub("0/1/2/3", "0/2", $i)}
				if (substr($i,1,8)=="0/1/2/3:" && substr($i,length($i)-11,7)=="0.0,0.0" && substr($i,length($i)-2,3)!="0.0")
					{sub("0/1/2/3", "0/3", $i)}
				}; print
			}
		else {print $0}
		}
	}
}' $input_file \
| awk \
'BEGIN{FS=OFS="\t"}
{
len=length($5) 
if ($5 !=  "." && len == 3)
	{for(i=1;i<=NF;i++)
		{if (substr($i,1,4)=="0/0:" && substr($i,length($i)-6,7)=="0.0,0.0")
			{sub("0.0,0.0", "0.0", $i)}
		if (substr($i,1,4)=="0/1:" && substr($i,length($i)-2,3)=="0.0")
			{sub(",0.0", "", $i)}
		if (substr($i,1,4)=="0/2:" && substr($i,length($i)-2,3)!="0.0")
			{sub("0.0,", "", $i)}
		}; print
	}
else
	{if ($5 !=  "." && len == 5)
		{for(i=1;i<=NF;i++)
			{if (substr($i,1,4)=="0/0:"  && substr($i,length($i)-7,8)==":0.0,0.0")
				{sub("0.0,0.0", "0.0", $i)}
			if (substr($i,1,4)=="0/0:" && substr($i,length($i)-10,11)=="0.0,0.0,0.0")
				{sub("0.0,0.0,0.0", "0.0", $i)}
			if (substr($i,1,4)=="0/1:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
				{sub(",0.0", "", $i)}
			if (substr($i,1,4)=="0/2:" && substr($i,length($i)-2,3)!="0.0")
				{sub("0.0,", "", $i)}
			if (substr($i,1,4)=="0/3:" && substr($i,length($i)-2,3)!="0.0" && substr($i,length($i)-8,4)==":0.0")
				{sub("0.0,", "", $i)}
			if (substr($i,1,4)=="0/2:"&& substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-6,3)!="0.0" && substr($i,length($i)-8,1)==":")
				{sub(",0.0", "", $i)}
			if (substr($i,1,4)=="0/1:" && substr($i,length($i)-6,7)=="0.0,0.0")
				{sub(",0.0,0.0", "", $i)}
			if (substr($i,1,4)=="0/2:" && substr($i,length($i)-2,3)=="0.0" && substr($i,length($i)-11,3)=="0.0" && substr($i,length($i)-6,3)!="0.0")
				{sub("0.0,", "", $i); sub(",0.0", "", $i)}
			if (substr($i,1,4)=="0/3:" && substr($i,length($i)-2,3)!="0.0")
				{sub("0.0,0.0,", "", $i)}
			}; print
		}
	else {print $0}
	}
}' \
| awk \
'BEGIN{FS=OFS="\t"}
{
ONE = 0; TWO = 0; THREE = 0
if ($0 ~ /\/1/)
	{ONE = 1}
if ($0 ~ /\/2/)
	{TWO = 1}
if ($0 ~ /\/3/)
	{THREE = 1}
if (ONE + TWO + THREE == 0)
	{sub($5,".", $5); print $0}
if (ONE + TWO + THREE >= 2)
	{next}
else
	{if (ONE == 1)
		{sub($5,substr($5,1,1),$5); print $0}
	else
		{if (TWO == 1)
			{sub($5,substr($5,3,1),$5); print $0}
		else
			{if (THREE == 1)
				{sub($5,substr($5,5,1),$5); print $0}
			}
		}
	}
}' \
| awk \
'BEGIN{FS=OFS="\t"} {gsub("/2","/1"); print $0}' \
| awk \
'BEGIN{FS=OFS="\t"} {gsub("/3","/1"); print $0}' > $output_file


