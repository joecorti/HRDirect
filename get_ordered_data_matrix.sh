set -e
set -u

if [ $# -ne 1 ]; then
	echo "Usage: $0 file" >&2
	echo "'file' with data-matrix, either HRDetect_fullPipeline.out or '*.data-matrix'" >&2
	exit 78
fi

##
## 6 fields required
##
##

#########################################################
function get_data_matrix {
# $1: HRDetect_fullPipeline.out file
	awk 'BEGIN{
			go_ahead=0
			N_header = 0
			N_Values = 0
			Expected_N_values = 7
		}
		{
			if(NR==1) {
					if(NF!=6) {
							go_ahead=1
						}
					$1=$1
					for(i=1;i<=NF;i++) {
							header[N_header] = $i;
							N_header++;
						}
				}
			if(NR==2) {
					for(j=1;j<=NF;j++) {
							values[N_Values] = $j;
							N_Values++;
						}
					if(N_Values==Expected_N_values) go_ahead=0
				}
			if(NR>2) {
					if(go_ahead==1) {
							if(NR%2==1) {
									$1=$1
									for(ii=1;ii<=NF;ii++) {
										header[N_header] = $ii;
										N_header++;
											}
								}
							else { # inizia dopo il nome del campione
								for(jj=2;jj<=NF;jj++) {
									values[N_Values] = $jj;
									N_Values++;
									}
								if(N_Values==Expected_N_values) go_ahead=0
								}
						}
					else exit
				}
		}
		END{
				header_string = ""
				for(idx in header) header_string = header_string " " header[idx]
				sub(/^ /, "", header_string)
				if(length(header_string)) {
					print header_string
					values_string = ""
					for(idx in values) values_string = values_string " " values[idx]
					sub(/^ /, "", values_string)
					print values_string
				}
			}' $1 | \
		_get_data_matrix -
}

function _get_data_matrix {
# stdin: output from 'get_raw_data_matrix_from_head'
	awk 'BEGIN{OFS="\t"}
		NR==1 {
			$1=$1
			# del.mh.prop SNV3 SV3      SV5 hrd SNV8
			for(i=1;i<=NF;i++) {
					if($i=="del.mh.prop") MH_i=i+1
					if($i=="SNV3") SNV3_i=i+1
					if($i=="SNV8") SNV8_i=i+1
					if($i=="SV3") SV3_i=i+1
					if($i=="SV5") SV5_i=i+1
					if($i=="hrd") hrd_i=i+1
				}

			print "", $(MH_i-1), $(SNV3_i-1), $(SV3_i-1), $(SV5_i-1), $(hrd_i-1), $(SNV8_i-1)
		}
		NR==2 {
			print $1, $MH_i, $SNV3_i, $SV3_i, $SV5_i, $hrd_i, $SNV8_i
		}'
}

#########################################################

get_data_matrix $1
