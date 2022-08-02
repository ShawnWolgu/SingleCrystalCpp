#include "singleX.h"

Matrix6d read_elastic(ifstream &is);
Matrix3d read_lattice(ifstream &is);
Matrix3d read_orientation(ifstream &is);
Matrix3d read_euler(ifstream &is);
void print_harden_law();
void add_slips(ifstream &is, vector<Slip> &slips, Matrix3d lattice_vecs);
Matrix3d vel_grad_flag_config(ifstream &load_file, Matrix3d &vel_grad_tensor);
Matrix3d load_matrix_input(ifstream &load_file);
void vel_grad_modify(char &flag_1, char &flag_2, int compoidx, Matrix3d &vel_grad_tensor);
void step_config(string in_str);
string get_next_line(ifstream &infile);
bool hasEnding (std::string const &fullString, std::string const &ending);

void set_config(){
    ifstream config;
    config.open(configure_path);
    if(config){
	cout << "Find configration file in cwd." << endl;
	string input_line;
	while (!config.eof()){
	    getline(config,input_line);
	    if (input_line[0] == '#'){
		if((input_line.find("Single") != input_line.npos) || (input_line.find("sx") != input_line.npos)){
		    sxfile_path = get_next_line(config);
		}
		if((input_line.find("Load") != input_line.npos) || (input_line.find("load") != input_line.npos)){
		    loadfile_path = get_next_line(config);
		}
	    }
	}
    }
}

Grain read_grain(){
    ifstream input_file(sxfile_path);
    if (!input_file) {
	cout << "Cannot find input_file!" << endl;
	exit(0);
    }
    string input_line;
    Matrix6d elastic_modulus;
    Matrix3d lattice_vecs, orientation;
    vector<Slip> slips;
    while (!input_file.eof())
    {
        getline(input_file, input_line);
        if(input_line[0] == '#'){
            if((input_line.find("Elastic") != input_line.npos) || (input_line.find("elastic") != input_line.npos)) {
                elastic_modulus = read_elastic(input_file);
                continue;
            }
            if((input_line.find("Harden") != input_line.npos) || (input_line.find("harden") != input_line.npos)) {
                getline(input_file,input_line);
                flag_harden = atoi(&input_line[0]);
                print_harden_law();
                continue;
            }
            if((input_line.find("Lattice") != input_line.npos) || (input_line.find("lattice") != input_line.npos)) {
                lattice_vecs = read_lattice(input_file);
                continue;
            }
            if((input_line.find("Orientation") != input_line.npos) || (input_line.find("orientation") != input_line.npos)) {
                orientation = read_orientation(input_file);
                continue;
            }
            if((input_line.find("Euler") != input_line.npos) || (input_line.find("euler") != input_line.npos)) {
                orientation = read_euler(input_file);
                continue;
            }
            if((input_line.find("Slip") != input_line.npos) || (input_line.find("slip") != input_line.npos)) {
                add_slips(input_file, slips, lattice_vecs); 
                continue;
            }
            continue;            
        }
        else{continue;}
    }
    input_file.close();
    for(auto &islip : slips) islip.cal_shear_modulus(elastic_modulus);
    return Grain(elastic_modulus, lattice_vecs, slips, orientation);
}

Matrix6d read_elastic(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[6] = {0,0,0,0,0,0};
    string temp_str;
    Matrix6d modulus;
    for(;row_num !=6; ++row_num){
        temp_idx = 0;
        getline(is, temp_str);
        stringstream stream(temp_str);
        while(!stream.eof()) stream >> temp[temp_idx++];
        modulus.row(row_num) << temp[0], temp[1], temp[2], temp[3], temp[4], temp[5];
    }
    cout << "Read Elastic Modulus:" << endl << modulus << endl;
    return modulus;
}

Matrix3d read_lattice(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[3] = {0,0,0};
    string temp_str;
    Matrix3d lattice;
    for(;row_num !=3; ++row_num){
        temp_idx = 0;
        getline(is, temp_str);
        stringstream stream(temp_str);
        while(!stream.eof()) stream >> temp[temp_idx++];
        lattice.row(row_num) << temp[0], temp[1], temp[2];
    }
    cout << "Read Lattice Vectors:" << endl << lattice << endl;
    return lattice;
}

Matrix3d read_orientation(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[3] = {0,0,0};
    string temp_str;
    Matrix3d orientation;
    for(;row_num !=3; ++row_num){
        temp_idx = 0;
        getline(is, temp_str);
        stringstream stream(temp_str);
        while(!stream.eof()) stream >> temp[temp_idx++];
        orientation.row(row_num) << temp[0], temp[1], temp[2];
    }
    cout << "Read Orientation Matrix:" << endl << orientation << endl;    
    return orientation;
}

Matrix3d read_euler(ifstream &is){
    int row_num = 0, temp_idx = 0;
    double temp[3] = {0,0,0};
    string temp_str;
    Vector3d euler_angle;
    Matrix3d orientation;
    temp_idx = 0;
    getline(is, temp_str);
    stringstream stream(temp_str);
    while(!stream.eof()) stream >> temp[temp_idx++];
    euler_angle << temp[0], temp[1], temp[2];
    cout << "Euler Angle:" << endl << euler_angle.transpose() << endl; 
    orientation = Euler_trans(euler_angle);
    cout << "Read Orientation Matrix (From Euler Angle):" << endl << orientation << endl; 
    return orientation;
}

void print_harden_law(){
    switch (flag_harden)
    {
    case 0:
        cout << "Hardening Law: Voce Hardening." << endl;
        break;
    case 1:
        cout << "Hardening Law: Dislocation Density Hardening." << endl;
        break;
    case 2:
        cout << "Hardening Law: Dislocation Velocity Model." << endl;
        break;
    default:
        cout << "Hardening Law: Voce Hardening." << endl;
        break;
    }
}

void add_slips(ifstream &is, vector<Slip> &slips, Matrix3d lattice_vecs){
    string input_line;
    vector<Vector6d> slip_infos;
    vector<double> harden_params;
    double temp;
    int slip_num = 0;

    getline(is, input_line);
    slip_num = atoi(&input_line[0]);
    cout << "Add " << slip_num << " Slips" << endl;
    for(int islip = 0; islip != slip_num; ++islip){
        getline(is, input_line);
        stringstream stream(input_line);
        Vector6d p_b;
        stream >> p_b(0) >> p_b(1) >> p_b(2) >> p_b(3) >> p_b(4) >> p_b(5);
        slip_infos.push_back(p_b);
    }
    getline(is, input_line);
    if(input_line[0] == '#')     getline(is, input_line);
    stringstream stream(input_line);
    while(stream >> temp) harden_params.push_back(temp);

    for(int islip = 0; islip != slip_num; ++islip){
        Slip temp_slip(slip_infos[islip], harden_params, lattice_vecs);
        cout << slip_infos[islip].transpose() << endl;
        for (auto i: harden_params) std::cout << i << ' ';
        cout << endl;
        slips.push_back(temp_slip);
    }
}

void read_load(Matrix3d &vel_grad_tensor, Matrix3d &vel_grad_flag, Matrix3d &stress_incr, Matrix3d &dstress_flag){
    cout << "Open Load File." << endl;
    ifstream load_file(loadfile_path,ifstream::in);
    if (!load_file) {
	cout << "Cannot find load_file!" << endl;
	exit(0);
    }
    string step_conf_string;
    getline(load_file, step_conf_string);
    // args: 1: timestep, 2: substep, 3: max_strain
    step_config(step_conf_string);
    cout << "Step configured." << endl;
    vel_grad_tensor = load_matrix_input(load_file);
    vel_grad_flag = vel_grad_flag_config(load_file,vel_grad_tensor);
    stress_incr = load_matrix_input(load_file);
    dstress_flag = load_matrix_input(load_file);
    if (vel_grad_flag.sum() + dstress_flag.sum() < 9) {throw "Cannot solve BC!";}
    cout << "Boundary conditions configured." << endl;
    cout << bc_modi_matrix << endl;
    load_file.close();
}

Matrix3d vel_grad_flag_config(ifstream &load_file, Matrix3d &vel_grad_tensor){
    if(load_file.eof()){
        throw "ERROR in load file!";
    }
    string matrix_str;
    Matrix3d temp_matrix = Matrix3d::Zero();
    char temp[9];
    char *temp_ar = temp;
    int line_count = 0;
    while(!load_file.eof() && line_count != 3){  
        getline(load_file,matrix_str);
        if(matrix_str == "" || matrix_str == "\r")  continue;
        stringstream stream(matrix_str);
        int temp_idx = 0;      
        for(; !stream.eof() && temp_idx != 3; temp_idx++) stream >> *(temp_ar++);
        if (temp_idx !=3) throw "ERROR in load file!";
        line_count++;
    }
    if (line_count != 3) throw "ERROR in load file!";
    if (temp[0] == 'w' || temp[4] == 'w' || temp[8] == 'w') throw "ERROR! w should not locate at the diagonal components!";
    vel_grad_modify(temp[5],temp[7], 3, vel_grad_tensor); 
    vel_grad_modify(temp[2],temp[6], 4, vel_grad_tensor);
    vel_grad_modify(temp[1],temp[3], 5, vel_grad_tensor);

    int temp_idx = 0;
    for (int i = 0;i!=3;++i){
	for (int j = 0;j != 3;++j){
	    unsigned char temp_n = temp[temp_idx++];
	    temp_matrix(i,j) = (float)(temp_n-'0');
	}
    }
    cout << vel_grad_tensor << endl;
    cout << temp_matrix << endl;
    return temp_matrix;
}

void vel_grad_modify(char &flag_1, char &flag_2, int compoidx, Matrix3d &vel_grad_tensor){
    // cases: 1, d, w
    int idx1 = 0, idx2 = 0;
    switch(compoidx){
	case 3: { idx1 = 1, idx2 = 2; break;}
	case 4: { idx1 = 0, idx2 = 2; break;}
	case 5: { idx1 = 0, idx2 = 1; break;}
    }
    if(flag_1 != flag_2 && (flag_1 != '0' && flag_2 != '0')){
	switch(flag_1){
	    case '1' :{switch(flag_2){
			case 'd': {vel_grad_tensor(idx2,idx1) = vel_grad_tensor(idx2,idx1) * 2 - vel_grad_tensor(idx1,idx2);break;}
			case 'w': {vel_grad_tensor(idx2,idx1) = vel_grad_tensor(idx1,idx2) - vel_grad_tensor(idx1,idx2) * 2;break;}
		      } break;}
	    case 'd' :{switch(flag_2){	
			case '1': {vel_grad_tensor(idx1,idx2) = vel_grad_tensor(idx1,idx2) * 2 - vel_grad_tensor(idx2,idx1);break;}
			case 'w': {float d = vel_grad_tensor(idx1,idx2), w = -vel_grad_tensor(idx2,idx1);
				vel_grad_tensor(idx1,idx2) = d + w; vel_grad_tensor(idx2,idx1) = d - w; break;}
		      } break;}
	    case 'w' :{switch(flag_2){	
			case '1': {vel_grad_tensor(idx1,idx2) = vel_grad_tensor(idx2,idx1) + 2 * vel_grad_tensor(idx1,idx2);break;}
			case 'd': {float d = vel_grad_tensor(idx2,idx1), w = vel_grad_tensor(idx1,idx2);
				vel_grad_tensor(idx1,idx2) = d + w; vel_grad_tensor(idx2,idx1) = d - w; break;}
		      } break;}
	}
	flag_1 = '1', flag_2 = '1'; 
    }
    Matrix<double,1,9> row_zero = Matrix<double,1,9>::Zero();
    if(flag_1 == flag_2){
	switch(flag_1){
	    case 'd' :{
		row_zero(0,compoidx) = 1; bc_modi_matrix.row(compoidx) = row_zero;
		row_zero(0,compoidx+3) = -1; bc_modi_matrix.row(compoidx+3) = row_zero; 
		flag_1 = '1', flag_2 = '0';
		} 
	    case 'w' :{
		row_zero(0,compoidx) = 1; bc_modi_matrix.row(compoidx+3) = row_zero;
		row_zero(0,compoidx+3) = 1; bc_modi_matrix.row(compoidx) = row_zero; 
		flag_1 = '1', flag_2 = '0';} 
	}}
    else{
	if(flag_1 == '0'){
	    Matrix<double,1,9> row_zero = Matrix<double,1,9>::Zero();
	    switch(flag_2){
		case 'd': {
		    row_zero(0,compoidx+3) = 1; bc_modi_matrix.row(compoidx) = row_zero;
		    row_zero(0,compoidx) = -1; bc_modi_matrix.row(compoidx+3) = -1 * row_zero;
		    break;} 
		case 'w': {
		    row_zero(0,compoidx+3) = -1; bc_modi_matrix.row(compoidx+3) = row_zero;
		    row_zero(0,compoidx) = -1; bc_modi_matrix.row(compoidx) = -1 * row_zero;
		    break;} 
		case '1': break;
	    }
	    flag_2 = '1';
	}
	else{
	    Matrix<double,1,9> row_zero = Matrix<double,1,9>::Zero();
	    switch(flag_1){
	    	case 'd' :{
		    row_zero(0,compoidx) = 1; bc_modi_matrix.row(compoidx) = row_zero;
		    row_zero(0,compoidx+3) = -1; bc_modi_matrix.row(compoidx+3) = row_zero; break;
		} 
		case 'w' :{
		    row_zero(0,compoidx) = 1; bc_modi_matrix.row(compoidx+3) = row_zero;
		    row_zero(0,compoidx+3) = 1; bc_modi_matrix.row(compoidx) = row_zero; break; 
		} 
		case '1': break;
	    }
	    flag_1 = '1';
	}
    }
}

Matrix3d load_matrix_input(ifstream &load_file){
    if(load_file.eof()){
        throw "ERROR in load file!";
    } 
    Matrix3d temp_matrix;
    string matrix_str;
    double temp[9]={0,0,0,0,0,0,0,0,0};
    double *temp_ar = temp;
    int line_count = 0;
    while(!load_file.eof() && line_count != 3){  
        getline(load_file,matrix_str);
        if(matrix_str == "" || matrix_str == "\r")  continue;
        stringstream stream(matrix_str);
        int temp_idx = 0;      
        for(; !stream.eof() && temp_idx != 3; temp_idx++) stream >> *(temp_ar++);
        if (temp_idx !=3) throw "ERROR in load file!";
        line_count++;
    }
    if (line_count != 3) throw "ERROR in load file!";
    temp_matrix << temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7], temp[8];
    cout << temp_matrix << endl;
    return temp_matrix;
}

void step_config(string in_str){
    stringstream stream(in_str);
    double temp[3]={0,0,0};
    int temp_idx = 0;
    while (!stream.eof()) stream >> temp[temp_idx++];
    timestep = 0.001, substep = 0.001, dtime = substep*timestep, max_strain = 0.2;

    if (temp[0] != 0) {
        timestep = temp[0];
    }

    if (temp[1] != 0) {
        substep = temp[1];
    }

    if (temp[2] != 0) {
        max_strain = temp[2];
    }
    dtime = substep*timestep;
}

bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

string get_next_line(ifstream &infile){
    string input_line;
    getline(infile,input_line);
    if(input_line == "" || input_line == "\r" || input_line[0] == '#') return get_next_line(infile);
    else{
    	if(hasEnding(input_line,"\r")){
	    input_line.erase(input_line.size() - 1);
	}
	return input_line;
    }
}

