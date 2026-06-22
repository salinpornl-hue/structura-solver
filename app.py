import streamlit as st
import pandas as pd
import subprocess
import os

st.title("Physical Results Analysis (Fortran + Streamlit)")

# ฟังก์ชันสำหรับคอมไพล์ Fortran (รันครั้งแรกครั้งเดียว)
@st.cache_resource
def compile_fortran():
    # ตรวจสอบชื่อไฟล์ให้ตรงกับบน GitHub ของคุณ
    fortran_file = "pointload2dynamics.f90" 
    if os.path.exists(fortran_file):
        try:
            # สั่งคอมไพล์โค้ดเป็นโปรแกรมชื่อ solver_app
            subprocess.run(["gfortran", "-o", "solver_app", fortran_file], check=True)
            return True
        except Exception as e:
            st.error(f"เกิดข้อผิดพลาดในการคอมไพล์: {e}")
            return False
    return False

# เริ่มการคอมไพล์
is_compiled = compile_fortran()

if is_compiled:
    st.success("คอมไพล์เครื่องจักรคำนวณ Fortran สำเร็จ!")
    
    # สร้างปุ่มให้กดคำนวณ
    if st.button("เริ่มคำนวณผลลัพธ์ (Run Engine)"):
        with st.spinner("Fortran กำลังคำนวณความเร็วสูง..."):
            try:
                # สั่งรันโปรแกรม Fortran เพื่อให้คายไฟล์ output_data.csv ออกมา
                subprocess.run(["./solver_app"], check=True)
                
                # ตรวจสอบว่าไฟล์สร้างเสร็จเรียบร้อยไหม
                if os.path.exists("output_data.csv"):
                    df = pd.read_csv("output_data.csv")
                    st.balloons() # เอฟเฟกต์แสดงความยินดี
                    
                    # แสดงตารางข้อมูล
                    st.subheader("ตารางข้อมูลผลลัพธ์")
                    st.dataframe(df)
                    
                    # วาดกราฟ
                    st.subheader("กราฟแสดงพฤติกรรม (Real Parts)")
                    st.line_chart(df.set_index("x")[["u_y_real", "u_x_real"]])
                else:
                    st.error("โปรแกรมรันผ่าน แต่ไม่พบไฟล์ output_data.csv ลองเช็คการเปิดไฟล์ในโค้ด Fortran")
            except Exception as e:
                st.error(f"รันโปรแกรมไม่สำเร็จ: {e}")
else:
    st.error("ไม่สามารถเตรียมระบบคำนวณได้ กรุณาตรวจสอบไฟล์ .f90 บน GitHub")
