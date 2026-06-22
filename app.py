import streamlit as st
import pandas as pd

st.title("Physical Results Analysis")

# อ่านไฟล์ CSV ที่ Fortran สร้างขึ้น
df = pd.read_csv("output_data.csv")

# แสดงตาราง
st.subheader("Data Table")
st.dataframe(df)

# แสดงกราฟเปรียบเทียบค่า
st.subheader("Visualization")
# ตัวอย่าง: Plot กราฟเส้นของ Real part ทั้งสองค่าเทียบกับ x
st.line_chart(df.set_index("x")[["u_y_real", "u_x_real"]])
