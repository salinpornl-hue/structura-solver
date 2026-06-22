import streamlit as st
import pandas as pd
import subprocess
import os

st.title("Physical Results Analysis (Fortran + Streamlit)")

# ฟังก์ชันสำหรับคอมไพล์ Fortran (เวอร์ชันอัปเกรดเพื่อรองรับ Matrix & ดักจับ Error)
@st.cache_resource
def compile_fortran():
    fortran_file = "pointload2dynamics.f90" 
    if os.path.exists(fortran_file):
        try:
            # เพิ่ม "-llapack" และ "-lblas" เพื่อให้ gfortran รู้จักคำสั่งคำนวณ Matrix
            # และใช้ capture_output เพื่อดึงข้อความ Error มาแสดงผล
            result = subprocess.run(
                ["gfortran", "-o", "solver_app", fortran_file, "-llapack", "-lblas"], 
                capture_output=True, 
                text=True
            )
            
            # ตรวจสอบว่าคอมไพล์สำเร็จไหม (ถ้า returncode ไม่เท่ากับ 0 แปลว่าพัง)
            if result.returncode != 0:
                st.error("❌ เกิดข้อผิดพลาดในโค้ด Fortran (Compilation Failed):")
                st.code(result.stderr, language="bash") # โชว์ Error ให้เห็นว่าบรรทัดไหนพัง
                return False
                
            return True
        except Exception as e:
            st.error(f"เกิดข้อผิดพลาดจากระบบ Python ในการสั่งคอมไพล์: {e}")
            return False
    else:
        st.error(f"ไม่พบไฟล์: {fortran_file} กรุณาตรวจสอบชื่อไฟล์บน GitHub")
        return False

# เริ่มการคอมไพล์
is_compiled = compile_fortran()

if is_compiled:
    st.success("🤖 คอมไพล์เครื่องจักรคำนวณ Fortran สำเร็จ!")
    
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
                    st.subheader("📊 ตารางข้อมูลผลลัพธ์")
                    st.dataframe(df)
                    
                    # วาดกราฟ
                    st.subheader("📈 กราฟแสดงพฤติกรรม (Real Parts)")
                    st.line_chart(df.set_index("x")[["u_y_real", "u_x_real"]])
                else:
                    st.error("⚠️ โปรแกรมรันผ่าน แต่ไม่พบไฟล์ output_data.csv ลองตรวจสอบตำแหน่งเขียนไฟล์ในโค้ด Fortran")
            except Exception as e:
                st.error(f"รันโปรแกรมไม่สำเร็จ: {e}")
else:
    st.error("🚫 ไม่สามารถเตรียมระบบคำนวณได้ กรุณาแก้ไขโค้ด Fortran ตามคำแนะนำด้านบน")
