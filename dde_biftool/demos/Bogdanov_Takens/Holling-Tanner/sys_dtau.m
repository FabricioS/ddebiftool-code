function dtau = sys_dtau(delay_nr,xx,par,nx,np)

dtau = [];

if length(nx) == 1 && isempty(np)
	switch nx
	case 0
		switch delay_nr
		case 1
			dtau = [0, 0];
		end
	case 1
		switch delay_nr
		case 1
			dtau = [0, 0];
		end
	end
elseif isempty(nx) && length(np) == 1
	switch np
	case 1
		switch delay_nr
		case 1
			dtau = 0;
		end
	case 2
		switch delay_nr
		case 1
			dtau = 0;
		end
	case 3
		switch delay_nr
		case 1
			dtau = 0;
		end
	case 4
		switch delay_nr
		case 1
			dtau = 0;
		end
	case 5
		switch delay_nr
		case 1
			dtau = 0;
		end
	case 6
		switch delay_nr
		case 1
			dtau = 0;
		end
	end
elseif length(nx) == 1 && length(np) == 1
	switch nx
	case 0
		switch np
		case 1
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 2
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 3
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 4
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 5
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 6
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		end
	case 1
		switch np
		case 1
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 2
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 3
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 4
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 5
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		case 6
			switch delay_nr
			case 1
				dtau = [0, 0];
			end
		end
	end
elseif length(nx) == 2 && isempty(np)
	nx1 = nx(1); nx2 = nx(2);
	switch nx1
	case 0
		switch nx2
		case 0
			switch delay_nr
			case 1
				dtau = [0, 0; 0, 0];
			end
		case 1
			switch delay_nr
			case 1
				dtau = [0, 0; 0, 0];
			end
		end
	case 1
		switch nx2
		case 0
			switch delay_nr
			case 1
				dtau = [0, 0; 0, 0];
			end
		case 1
			switch delay_nr
			case 1
				dtau = [0, 0; 0, 0];
			end
		end
	end
end
if isempty(dtau)
	display([delay_nr nx np]);
	error('SYS_DTAU: requested derivative could not be computed!');
end