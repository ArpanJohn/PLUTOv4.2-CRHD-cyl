#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                CYLINDRICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     8

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO
#define  CR_FLUID                NC_PdV_TOTENG
#define  CR_DIFFUSION            SUPER_TIME_STEPPING

/* -- user-defined parameters (labels) -- */

#define  E_SN                    0
#define  M_SN                    1
#define  RHO_AMB                 2
#define  TEMP_AMB                3
#define  mu_AMB                  4
#define  Pcr_Shock               5
#define  OPT                     6
#define  Rinj                    7

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH             CONST_pc
#define  UNIT_VELOCITY           (1e5)
#define  UNIT_DENSITY            CONST_mH
#define  unitPRS                 (UNIT_DENSITY*pow(UNIT_VELOCITY,2.0))
#define  unitMASS                (UNIT_DENSITY*pow(UNIT_LENGTH,3.0))
#define  unitTIME                (UNIT_LENGTH/UNIT_VELOCITY)
#define  unitENERGY              (unitMASS*pow(UNIT_LENGTH,2.0)*pow(unitTIME,-2.0))

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             YES
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   DEFAULT
#define  CT_EMF_AVERAGE            UCT_HLL
#define  CT_EN_CORRECTION          NO
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO
