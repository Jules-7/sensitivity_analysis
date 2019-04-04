import os
from general.constants import ft, knot, ton


class AircraftPerformanceModel(object):

    def __init__(self, actype):
        self.actype = actype
        self.path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../BADA_3_12")) + "/"
        self.extract_parameters()
        return

    def extract_parameters(self):

        # ----------- BADA User Manual rev 3-12. p.33-36 ------
        # ------------ GLOBAL AIRCRAFT PARAMETERS -------------

        # ------------ MAX ACCELERATION (5.2) -----------------------
        self.acceler_long_max = 2.0 * ft  # [m/s**2] max longitudinal acceleration for TAS (2 ft/s**2)
        self.acceler_norm_max = 5.0 * ft  # [m/s**2] max normal acceleration for ROCD (5 ft/s**2)

        # ------------ BANK ANGLES (5.3) ----------------------------
        self.phi_nom_to_ld = 15  # [deg] Nominal bank angles for civil flight during TO and LD
        self.phi_nom = 35  # [deg] Nominal bank angles for civil flight during all other phases

        self.phi_max_to_ld = 25  # [deg] Maximum bank angles for civil flight during TO and LD
        self.phi_max_hold = 35  # [deg] Maximum bank angles for civil flight during HOLD
        self.phi_max = 45  # [deg] Maximum bank angles for civil flight during all other phases

        # ----------- THRUST FACTORS (5.5) -------------------------------------
        self.C_th_to = 1.2  # [-] Take-off thrust coefficient Maximum
        self.C_t_cr = 0.95  # [-] Maximum cruise thrust coefficient p.34, used in engine_thrust module

        # ------------ CONFIGURATION ALTITUDE THRESHOLD (5.6) ---------------------------
        self.H_max_to = 121.92  # [m] Maximum altitude threshold for take-off - (400 ft)
        self.H_max_ic = 609.6   # [m] Maximum altitude threshold for initial climb - (2,000 ft)
        self.H_max_ap = 2438.4  # [m] Maximum altitude threshold for approach - (8,000 ft)
        self.H_max_ld = 914.4   # [m] Maximum altitude threshold for landing - (3,000 ft)

        # ---------- MIN SPEED COEFFICIENTS (5.7) -----------------------------
        # to be used in in formulas (4.1-*) and (4.3-*)
        self.C_v_min = 1.3  # [-] Minimum speed coefficient (all other phases)
        self.C_v_min_to = 1.2  # [-] Minimum speed coefficient for take-off

        # ----------- SPEED SCHEDULES (5.8) ------------------------------------
        # speeds converted into m/s
        self.V_d_cl_1 = 2.57222   # [m/s] Climb speed increment below FL15 (jet) 5 kts
        self.V_d_cl_2 = 5.14444   # [m/s] Climb speed increment below FL30 (jet) 10 kts
        self.V_d_cl_3 = 15.43332  # [m/s] Climb speed increment below FL40 (jet) 30 kts
        self.V_d_cl_4 = 30.86664  # [m/s] Climb speed increment below FL50 (jet) 60 kts
        self.V_d_cl_5 = 41.15552  # [m/s] Climb speed increment below FL60 (jet) 80 kts

        self.V_d_des_1 = 2.57222   # [m/s] DescentPrediction speed increment below FL10 (jet/turboprop) 5 kts
        self.V_d_des_2 = 5.14444   # [m/s] DescentPrediction speed increment below FL15 (jet/turboprop) 10 kts
        self.V_d_des_3 = 10.28888  # [m/s] DescentPrediction speed increment below FL20 (jet/turboprop) 20 kts
        self.V_d_des_4 = 25.7222   # [m/s] DescentPrediction speed increment below FL30 (jet/turboprop) 50 kts

        # ----------- REDUCED POWER COEFFICIENT (5.11) -------------------------------------
        self.C_red_ac = 0.15  # Maximum reduction in power for jets

        # ------------ END GLOBAL AIRCRAFT PARAMETERS -------------

        # ------------- READ OPF FILE -----------------------------
        file_name = self.path + self.actype + "__.OPF"

        with open(file_name) as opf:
            opf_data = opf.readlines()
            opf_data = [each.split() for each in opf_data]
            # opf_data = [' '.join(filter(None, each)) for each in opf_data]
            # opf_data = [each.split(" ") for each in opf_data]

            # ----------- BADA User Manual rev 3-12. p.27-28 ------
            # ---------- AIRCRAFT TYPE -----------------------
            self.n_eng = int(opf_data[13][2])  # number of engines
            self.engine_type = opf_data[13][4]  # Jet, Turboprop or Piston
            self.wake_category = opf_data[13][5]  # J, H, M, L
            # ---------- END AIRCRAFT TYPE ------------------

            # ----------- MASS ------------------------------
            self.mass_ref = float(opf_data[18][1]) * ton   # [kg]
            self.mass_min = float(opf_data[18][2]) * ton   # [kg]
            self.mass_max = float(opf_data[18][3]) * ton   # [kg]
            self.mass_pyld = float(opf_data[18][4]) * ton  # [kg] max payload mass
            self.G_w = float(opf_data[18][5]) * ft  # [ft/kg] weight gradient on max altitude
            # ---------- END MASS --------------------------

            # ----------- FLIGHT ENVELOPE -----------------
            self.V_mo = float(opf_data[21][1]) * knot  # [m/s] CAS max operating speed
            self.M_mo = float(opf_data[21][2])  # max operating Mach number
            self.h_mo = float(opf_data[21][3]) * ft  # [m] max operating altitude
            self.h_max = float(opf_data[21][4]) * ft  # [m] max altitude at MTOW and ISA
            self.G_t = float(opf_data[21][5]) * ft  # [ft/K] temp gradient on max altitude
            # ------------ END FLIGHT ENVELOPE ---------------

            # ------------ AERODYNAMICS ----------------------
            self.S = float(opf_data[25][2])  # [m**2] reference wing surface area
            self.C_lbo = float(opf_data[25][3])  # Buffet onset lift coef. (jet only)
            self.k = float(opf_data[25][4])  # buffering gradient (jet only)

            self.V_stall_cr_ref = float(opf_data[28][4]) * knot  # [m/s] CAS stall speed TO
            self.C_D0_cr = float(opf_data[28][5])  # parasitic drag coefficient TO
            self.C_D2_cr = float(opf_data[28][6])  # induced drag coefficient TO

            self.V_stall_ic_ref = float(opf_data[29][4]) * knot  # [m/s] CAS stall speed IC
            self.C_D0_ic = float(opf_data[29][5])  # parasitic drag coefficient IC
            self.C_D2_ic = float(opf_data[29][6])  # induced drag coefficient IC

            self.V_stall_to_ref = float(opf_data[30][4]) * knot  # [m/s] CAS stall speed CR
            self.C_D0_to = float(opf_data[30][5])  # parasitic drag coefficient (cruise)
            self.C_D2_to = float(opf_data[30][6])  # induced drag coefficient (cruise)

            self.V_stall_ap_ref = float(opf_data[31][4]) * knot  # [m/s] CAS stall speed AP
            self.C_D0_ap = float(opf_data[31][5])  # parasitic drag coefficient (approach)
            self.C_D2_ap = float(opf_data[31][6])  # induced drag coefficient (approach)

            self.V_stall_ld_ref = float(opf_data[32][4]) * knot  # [m/s] CAS stall speed LD
            self.C_D0_ld = float(opf_data[32][5])  # parasitic drag coefficient (landing)
            self.C_D2_ld = float(opf_data[32][6])  # induced drag coefficient (landing)

            self.C_D0_delta_ldg = float(opf_data[38][3])  # parasite drag coef. (landing gear)
            # ------------ END AERODYNAMICS ----------------------

            # -------------- ENGINE THRUST --------------------
            self.C_tc_1 = float(opf_data[44][1])  # [N] 1st max. climb thrust coefficient
            self.C_tc_2 = float(opf_data[44][2]) * ft  # [m] 2nd max climb thrust coefficient
            self.C_tc_3 = float(opf_data[44][3]) * 1./ft**2  # [1/m**2] 3rd max. climb thrust coefficient
            self.C_tc_4 = float(opf_data[44][4])  # [K] 1st thrust temperature coefficient
            self.C_tc_5 = float(opf_data[44][5])  # [1/K] 2nd thrust temperature coefficient

            self.C_t_des_low = float(opf_data[46][1])  # low altitude descent thrust coefficient
            self.C_t_des_high = float(opf_data[46][2])  # high altitude descent thrust coefficient
            self.H_p_des = float(opf_data[46][3]) * ft  # [m] transition altitude for calculation of descent thrust
            self.C_t_des_app = float(opf_data[46][4])  # approach thrust coefficient
            self.C_t_des_ld = float(opf_data[46][5])  # landing thrust coefficient

            self.V_des_ref = float(opf_data[48][1]) * knot  # [m/s] reference descent speed (CAS)

            self.M_des_ref = float(opf_data[48][2])  # reference descent Mach number
            # -------------- END ENGINE THRUST ----------------

            # -------------- FUEL CONSUMPTION -----------------
            self.C_f_1 = float(opf_data[51][1])  # 1st thrust specific fuel consumption coefficient
            self.C_f_2 = float(opf_data[51][2])  # 2nd thrust specific fuel consumption coefficient

            self.C_f_3 = float(opf_data[53][1])  # 1st descent fuel flow coefficient
            self.C_f_4 = float(opf_data[53][2])  # 2nd descent fuel flow coefficient

            self.C_f_cr = float(opf_data[55][1])  # cruise fuel flow correction coefficient
            # -------------- END FUEL CONSUMPTION -----------------

        # ------------- END READ OPF FILE -----------------------------

        # ------------- READ APF FILE -----------------------------
        file_name = self.path + self.actype + "__.APF"

        with open(file_name) as apf:
            apf_data = apf.readlines()
            apf_data = [each.split() for each in apf_data]

            if apf_data[21][2] == "AV":
                line_number = 21
            else:  # if no data for average weight - take low
                line_number = 20
            # ------------------ CLIMB --------------------------------
            # ---------------- BADA User Manual rev. 3-12, p.30
            # parameters to characterise the climb phase
            self.V_cl_1 = float(apf_data[line_number][3]) * knot  # [m/s] standard climb CAS between 1,500/6,000 and 10,000 ft
            self.V_cl_2 = float(apf_data[line_number][4]) * knot  # [m/s] standard climb CAS [knots] between 10,000 ft and Mach transition altitude
            self.M_cl = float(apf_data[line_number][5]) / 100.  # standard climb Mach number above Mach transition altitude

            # -------------- CRUISE ---------------------------------------------------------
            # -------------- BADA User Manual rev. 3-12, p.31
            # parameters to characterise the cruise phase
            self.V_cr_1 = float(apf_data[line_number][6]) * knot  # [m/s] standard cruise CAS between 3,000 and 10,000 ft
            self.V_cr_2 = float(apf_data[line_number][7]) * knot  # [m/s] standard cruise CAS between 10,000 ft and Mach transition altitude Mcr
            self.M_cr = float(apf_data[line_number][8]) / 100.  # standard cruise Mach number above Mach transition altitude

            # ---------- DESCENT ---------------------------------------------------------
            # -------------- BADA User Manual rev. 3-12, p.32
            # parameters to characterise the descent phase
            self.M_des = float(apf_data[line_number][9]) / 100.  # standard descent Mach number above Mach transition altitude
            self.V_des_1 = float(apf_data[line_number][10]) * knot  # [m/s] standard descent CAS [knots] between 3,000/6,000 and 10,000 ft
            self.V_des_2 = float(apf_data[line_number][11]) * knot  # [m/s] standard descent CAS [knots] between 10,000 ft and Mach transition altitude
        # ------------- END READ APF FILE -----------------------------
        return


