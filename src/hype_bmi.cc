#include <cstring>
#include "hype_bmi.h"


extern "C" int init_hype(char*, int*);
extern "C" int update_hype();
extern "C" int finalize_hype();
extern "C" int get_num_basins();
//extern "C" int get_hype_time(double*);
extern "C" int get_num_ovars();
extern "C" void get_ovar(char*, int*);

HypeBmi::HypeBmi(){}

HypeBmi::~HypeBmi(){}

int HypeBmi::initialize(const char *config_file)
{
    int mem = 0;
    init_hype(const_cast<char*>(config_file), &mem);
    return BMI_SUCCESS;
}

int HypeBmi::update()
{
    int status = update_hype();
    if(status == 0)
    {
        return BMI_SUCCESS;
    }
    return BMI_FAILURE;
}

int HypeBmi::update_until(double time)
{
    double t;
    int status = this->get_current_time(&t);
    while(status == BMI_SUCCESS and t < time)
    {
        this->update();
        status = this->get_current_time(&t);
    }
    return status;
}

int HypeBmi::update_frac(double time)
{
    return BMI_FAILURE;
}

int HypeBmi::run_model()
{
    return BMI_FAILURE;
}

int HypeBmi::finalize()
{
    finalize_hype();
    return BMI_SUCCESS;
}

int HypeBmi::get_component_name(char* name) const
{
    strcpy(name, "hype-5.6.1\0");
    return BMI_SUCCESS;
}

int HypeBmi::get_input_var_name_count(int* count) const
{
    *count = 0;
    return BMI_SUCCESS;
}

int HypeBmi::get_output_var_name_count(int* count) const
{
    *count = get_num_ovars();
    return BMI_SUCCESS;
}

int HypeBmi::get_input_var_names(char **names) const
{
    return BMI_SUCCESS;
}

int HypeBmi::get_output_var_names(char **names) const
{
   int m = get_num_ovars();
   int k = BMI_MAX_VAR_NAME;
   char* store = (char*)malloc(m * k * sizeof(char));
   char** x = (char**)malloc(m * sizeof(char*));
   for(int i = 0; i < m; ++i)
   {
        x[i] = &(store[i * k]);
   }
   get_ovar(store, &k);
   for(int i = 0; i < m; ++i)
   {
        strcpy(names[i], x[i]);
   }
   free(x);
   free(store);
   return BMI_SUCCESS;
}

int HypeBmi::get_var_type(const char* name, char* otype) const
{
   int m = get_num_ovars();
   int k = BMI_MAX_VAR_NAME;
   char* store = (char*)malloc(m * k * sizeof(char));
   // TODO: Add input variables whenever they are available
   get_ovar(store, &k);
   for(int i = 0; i < m; ++i)
   {
        if(strcmp(name, &(store[i * k])) == 0)
        {
            strcpy(otype, "double\n");
            return BMI_SUCCESS;
        }
   }
   free(store);
   return BMI_FAILURE;
}

int HypeBmi::get_var_grid(const char* name, int* gtype) const
{
   int m = get_num_ovars();
   int k = BMI_MAX_VAR_NAME;
   char* store = (char*)malloc(m * k * sizeof(char));
   // TODO: Add input variables whenever they are available
   get_ovar(store, &k);
   for(int i = 0; i < m; ++i)
   {
        if(strcmp(name, &(store[i * k])) == 0)
        {
            *gtype = 1;
            return BMI_SUCCESS;
        }
   }
   free(store);
   return BMI_FAILURE;
}


int HypeBmi::get_var_itemsize(const char* name, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_var_units(const char* name, char* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_var_nbytes(const char* name, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_current_time(double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_start_time(double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_end_time(double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_time_units(char* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_time_step(double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_value(const char* name, void* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_value_ptr(const char* name, void** dest)
{
   return BMI_FAILURE;
}

int HypeBmi::get_value_at_indices(const char* name, void* dest, const int* pts, int numpts) const
{
   return BMI_FAILURE;
}

int HypeBmi::set_value(const char* name, const void* src)
{
   return BMI_FAILURE;
}

int HypeBmi::set_value_ptr(const char* name, void** src)
{
   return BMI_FAILURE;
}

int HypeBmi::set_value_at_indices(const char* name, const int* pts, int numpts, const void* src)
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_size(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_rank(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_type(int id, char* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_shape(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_spacing(int id, double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_origin(int id, double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_x(int id, double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_y(int id, double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_z(int id, double* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_cell_count(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_point_count(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_vertex_count(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_connectivity(int id, int* dest) const
{
   return BMI_FAILURE;
}

int HypeBmi::get_grid_offset(int id, int* dest) const
{
   return BMI_FAILURE;
}
