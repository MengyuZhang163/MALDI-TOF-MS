import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
import zipfile
import io
import tempfile
import subprocess
import os
import shutil

# 页面配置
st.set_page_config(
    page_title="MALDI-TOF MS 数据处理平台",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 自定义CSS样式
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: 700;
        color: #1f77b4;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #555;
        margin-bottom: 2rem;
    }
    .stAlert {
        border-radius: 10px;
    }
    .upload-box {
        border: 2px dashed #1f77b4;
        border-radius: 10px;
        padding: 2rem;
        text-align: center;
        background-color: #f0f8ff;
    }
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 8px;
        padding: 1rem;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
</style>
""", unsafe_allow_html=True)

# 初始化session state
if 'processed_data' not in st.session_state:
    st.session_state.processed_data = None
if 'processing_params' not in st.session_state:
    st.session_state.processing_params = {}

def create_r_script(temp_dir, params):
    """生成自适应参数的R处理脚本"""
    r_script = f"""
library('MALDIquant')
library('MALDIquantForeign')
library('readxl')

# 自适应参数估计函数
estimate_halfWindowSize <- function(spectra, test_sizes = c(10, 20, 50, 90, 150)) {{
  peak_counts <- sapply(test_sizes, function(hw) {{
    test_peaks <- detectPeaks(spectra[[1]], method = "MAD", 
                               halfWindowSize = hw, SNR = 2)
    length(test_peaks@mass)
  }})
  optimal_idx <- which.min(abs(peak_counts - median(peak_counts)))
  return(test_sizes[optimal_idx])
}}

estimate_SNR <- function(spectra, quantile_threshold = 0.95) {{
  all_intensities <- unlist(lapply(spectra, function(s) s@intensity))
  noise_level <- quantile(all_intensities, 0.5)
  signal_level <- quantile(all_intensities, quantile_threshold)
  snr <- signal_level / noise_level
  return(max(2, min(5, round(snr * 0.1))))
}}

estimate_tolerance <- function(spectra) {{
  mz_range <- range(spectra[[1]]@mass)
  resolution <- length(spectra[[1]]@mass) / diff(mz_range)
  rel_tolerance <- 1 / resolution * 10
  return(max(0.002, min(0.01, rel_tolerance)))
}}

# 路径设置
train_path <- '{temp_dir}/train/'
valid_path <- '{temp_dir}/valid/'
excel_file <- list.files(train_path, pattern = '\\\\.xlsx$', full.names = TRUE)[1]

# 读取样本信息
samples <- read_excel(excel_file)

# 导入训练集
training_spectra <- importTxt(train_path)
cat(sprintf("导入了 %d 个训练集光谱\\n", length(training_spectra)))

# 参数估计
{'optimal_hw <- ' + str(params['halfWindowSize']) if params['auto_params'] == False else 'optimal_hw <- estimate_halfWindowSize(training_spectra)'}
{'optimal_snr <- ' + str(params['SNR']) if params['auto_params'] == False else 'optimal_snr <- estimate_SNR(training_spectra)'}
{'optimal_tol <- ' + str(params['tolerance']) if params['auto_params'] == False else 'optimal_tol <- estimate_tolerance(training_spectra)'}

cat(sprintf("参数设置: halfWindowSize=%d, SNR=%.1f, tolerance=%.4f\\n", 
            optimal_hw, optimal_snr, optimal_tol))

# 保存参数到文件
params_df <- data.frame(
  parameter = c('halfWindowSize', 'SNR', 'tolerance'),
  value = c(optimal_hw, optimal_snr, optimal_tol)
)
write.csv(params_df, '{temp_dir}/parameters.csv', row.names = FALSE)

# 预处理训练集
training_spectra <- transformIntensity(training_spectra, method = "sqrt")
training_spectra <- smoothIntensity(training_spectra, method = "SavitzkyGolay", 
                                     halfWindowSize = optimal_hw)
training_spectra <- removeBaseline(training_spectra, method = "SNIP", iterations = 100)
training_spectra <- calibrateIntensity(training_spectra, method = "TIC")

# 分配标签
train_labels <- samples$group[match(
  sapply(training_spectra, function(s) basename(s@metaData$file)),
  samples$file
)]

# 计算平均谱
avgSpectra <- averageMassSpectra(training_spectra, labels = train_labels)
avgSpectra <- alignSpectra(avgSpectra, halfWindowSize = optimal_hw,
                           SNR = optimal_snr, tolerance = optimal_tol,
                           warpingMethod = "lowess")

# 处理验证集（如果存在）
if (dir.exists(valid_path) && length(list.files(valid_path, pattern = '\\\\.txt$')) > 0) {{
  validation_spectra <- importTxt(valid_path)
  cat(sprintf("导入了 %d 个验证集光谱\\n", length(validation_spectra)))
  
  validation_spectra <- transformIntensity(validation_spectra, method = "sqrt")
  validation_spectra <- smoothIntensity(validation_spectra, method = "SavitzkyGolay", 
                                         halfWindowSize = optimal_hw)
  validation_spectra <- removeBaseline(validation_spectra, method = "SNIP", iterations = 100)
  validation_spectra <- calibrateIntensity(validation_spectra, method = "TIC")
  
  combinedSpectra <- c(avgSpectra, validation_spectra)
}} else {{
  combinedSpectra <- avgSpectra
  validation_spectra <- list()
}}

# 对齐、检峰、分箱
alignedCombined <- alignSpectra(combinedSpectra, halfWindowSize = optimal_hw,
                                SNR = optimal_snr, tolerance = optimal_tol,
                                warpingMethod = "lowess")

detectedPeaksCombined <- detectPeaks(alignedCombined, method = "MAD",
                                     halfWindowSize = optimal_hw, SNR = optimal_snr)

binnedPeaksCombined <- binPeaks(detectedPeaksCombined, tolerance = 2)

mat <- intensityMatrix(binnedPeaksCombined, alignedCombined)

# 处理列名
bin_centers_highprec <- as.numeric(colnames(mat))
bin_centers_integer <- round(bin_centers_highprec)
colnames(mat) <- paste0("mz_", bin_centers_integer)

# 处理行名
n_train_groups <- length(avgSpectra)
n_val_spectra <- length(validation_spectra)
train_group_names <- unique(train_labels)

if (n_val_spectra > 0) {{
  valid_names <- sapply(validation_spectra, function(s) basename(s@metaData$file))
  rownames(mat) <- c(train_group_names, valid_names)
}} else {{
  rownames(mat) <- train_group_names
}}

# 保存结果
intensity_train <- mat[1:n_train_groups, , drop = FALSE]
write.csv(intensity_train, '{temp_dir}/peak_intensity_train.csv', row.names = TRUE)

if (n_val_spectra > 0) {{
  intensity_validation <- mat[(n_train_groups + 1):(n_train_groups + n_val_spectra), , drop = FALSE]
  write.csv(intensity_validation, '{temp_dir}/peak_intensity_validation.csv', row.names = TRUE)
}}

# 保存第一个平均谱用于可视化
first_avg <- avgSpectra[[1]]
viz_data <- data.frame(
  mz = first_avg@mass,
  intensity = first_avg@intensity
)
write.csv(viz_data, '{temp_dir}/spectrum_viz.csv', row.names = FALSE)

cat("处理完成!\\n")
"""
    
    script_path = Path(temp_dir) / "process.R"
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(r_script)
    
    return script_path

def run_r_script(script_path):
    """执行R脚本"""
    try:
        result = subprocess.run(
            ['Rscript', str(script_path)],
            capture_output=True,
            text=True,
            timeout=300  # 5分钟超时
        )
        return result.stdout, result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        return "", "处理超时（超过5分钟）", 1
    except FileNotFoundError:
        return "", "未找到R环境，请确保已安装R和必要的包", 1

def plot_spectrum(df):
    """绘制质谱图"""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df['mz'],
        y=df['intensity'],
        mode='lines',
        name='强度',
        line=dict(color='#1f77b4', width=1),
        fill='tozeroy',
        fillcolor='rgba(31, 119, 180, 0.3)'
    ))
    
    fig.update_layout(
        title='平均质谱图',
        xaxis_title='m/z',
        yaxis_title='相对强度',
        hovermode='x unified',
        template='plotly_white',
        height=500
    )
    
    return fig

def plot_heatmap(df):
    """绘制强度热图"""
    # 取前50个最大峰
    top_cols = df.iloc[:, 1:].sum().nlargest(50).index
    data_subset = df[['行名'] + list(top_cols)]
    
    fig = px.imshow(
        data_subset.set_index('行名').T,
        aspect='auto',
        color_continuous_scale='Viridis',
        labels=dict(x="样本", y="m/z", color="强度")
    )
    
    fig.update_layout(
        title='峰强度热图（Top 50峰）',
        height=600
    )
    
    return fig

# ========================================
# 主应用界面
# ========================================

st.markdown('<div class="main-header">🔬 MALDI-TOF MS 数据处理平台</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">微生物质谱数据自动化预处理工具</div>', unsafe_allow_html=True)

# 侧边栏
with st.sidebar:
    st.header("📋 处理流程")
    st.markdown("""
    1️⃣ 上传训练集文件  
    2️⃣ 上传验证集文件（可选）  
    3️⃣ 配置处理参数  
    4️⃣ 开始处理  
    5️⃣ 查看结果并下载  
    """)
    
    st.divider()
    
    st.header("⚙️ 参数配置")
    
    auto_params = st.checkbox("自动参数估计", value=True, 
                              help="根据数据特征自动选择最佳参数")
    
    if not auto_params:
        st.subheader("手动参数设置")
        halfWindowSize = st.slider("半峰宽", 10, 200, 90, 10)
        SNR = st.slider("信噪比阈值", 1.0, 10.0, 2.0, 0.5)
        tolerance = st.slider("对齐容差", 0.001, 0.02, 0.008, 0.001, format="%.4f")
    else:
        halfWindowSize, SNR, tolerance = None, None, None
    
    st.session_state.processing_params = {
        'auto_params': auto_params,
        'halfWindowSize': halfWindowSize,
        'SNR': SNR,
        'tolerance': tolerance
    }
    
    st.divider()
    
    st.markdown("""
    ### 💡 使用提示
    - 训练集需包含TXT光谱文件和Excel标签文件
    - Excel必须有`file`和`group`两列
    - 验证集仅需TXT文件
    """)

# 主内容区
tab1, tab2, tab3 = st.tabs(["📁 数据上传", "▶️ 处理与结果", "📊 数据可视化"])

with tab1:
    st.header("数据上传")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("训练集文件")
        train_txt_files = st.file_uploader(
            "上传训练集TXT文件",
            type=['txt'],
            accept_multiple_files=True,
            key='train_txt'
        )
        
        train_excel = st.file_uploader(
            "上传标签Excel文件",
            type=['xlsx', 'xls'],
            key='train_excel',
            help="必须包含 'file' 和 'group' 两列"
        )
        
        if train_txt_files and train_excel:
            st.success(f"✅ 已上传 {len(train_txt_files)} 个TXT文件和1个Excel文件")
    
    with col2:
        st.subheader("验证集文件（可选）")
        valid_txt_files = st.file_uploader(
            "上传验证集TXT文件",
            type=['txt'],
            accept_multiple_files=True,
            key='valid_txt'
        )
        
        if valid_txt_files:
            st.success(f"✅ 已上传 {len(valid_txt_files)} 个验证集文件")
        else:
            st.info("💡 验证集为可选项")

with tab2:
    st.header("数据处理")
    
    if st.button("🚀 开始处理", type="primary", use_container_width=True):
        if not train_txt_files or not train_excel:
            st.error("❌ 请先上传训练集文件！")
        else:
            with st.spinner("正在处理数据..."):
                # 创建临时目录
                temp_dir = tempfile.mkdtemp()
                train_dir = Path(temp_dir) / "train"
                valid_dir = Path(temp_dir) / "valid"
                train_dir.mkdir()
                valid_dir.mkdir()
                
                try:
                    # 保存训练集文件
                    for txt_file in train_txt_files:
                        with open(train_dir / txt_file.name, 'wb') as f:
                            f.write(txt_file.read())
                    
                    with open(train_dir / train_excel.name, 'wb') as f:
                        f.write(train_excel.read())
                    
                    # 保存验证集文件
                    if valid_txt_files:
                        for txt_file in valid_txt_files:
                            with open(valid_dir / txt_file.name, 'wb') as f:
                                f.write(txt_file.read())
                    
                    # 生成并执行R脚本
                    script_path = create_r_script(temp_dir, st.session_state.processing_params)
                    
                    stdout, stderr, returncode = run_r_script(script_path)
                    
                    if returncode == 0:
                        st.success("✅ 处理完成！")
                        
                        # 读取结果
                        results = {}
                        
                        if (Path(temp_dir) / "peak_intensity_train.csv").exists():
                            results['train'] = pd.read_csv(Path(temp_dir) / "peak_intensity_train.csv")
                            results['train'].rename(columns={results['train'].columns[0]: '行名'}, inplace=True)
                        
                        if (Path(temp_dir) / "peak_intensity_validation.csv").exists():
                            results['validation'] = pd.read_csv(Path(temp_dir) / "peak_intensity_validation.csv")
                            results['validation'].rename(columns={results['validation'].columns[0]: '行名'}, inplace=True)
                        
                        if (Path(temp_dir) / "spectrum_viz.csv").exists():
                            results['spectrum'] = pd.read_csv(Path(temp_dir) / "spectrum_viz.csv")
                        
                        if (Path(temp_dir) / "parameters.csv").exists():
                            results['params'] = pd.read_csv(Path(temp_dir) / "parameters.csv")
                        
                        st.session_state.processed_data = results
                        
                        # 显示处理信息
                        st.subheader("处理摘要")
                        
                        col1, col2, col3 = st.columns(3)
                        
                        with col1:
                            st.metric("训练集样本数", len(results['train']))
                        
                        with col2:
                            if 'validation' in results:
                                st.metric("验证集样本数", len(results['validation']))
                            else:
                                st.metric("验证集样本数", "N/A")
                        
                        with col3:
                            st.metric("检测到的峰数", len(results['train'].columns) - 1)
                        
                        # 显示参数
                        if 'params' in results:
                            st.subheader("使用的处理参数")
                            st.dataframe(results['params'], use_container_width=True)
                        
                        # 显示日志
                        with st.expander("查看处理日志"):
                            st.code(stdout, language='text')
                    
                    else:
                        st.error(f"❌ 处理失败！\n\n{stderr}")
                        st.code(stdout, language='text')
                
                except Exception as e:
                    st.error(f"❌ 发生错误: {str(e)}")
                
                finally:
                    # 清理临时文件
                    shutil.rmtree(temp_dir, ignore_errors=True)
    
    st.divider()
    
    # 结果下载
    if st.session_state.processed_data:
        st.subheader("📥 下载处理结果")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if 'train' in st.session_state.processed_data:
                csv_train = st.session_state.processed_data['train'].to_csv(index=False)
                st.download_button(
                    label="下载训练集结果 (CSV)",
                    data=csv_train,
                    file_name="peak_intensity_train.csv",
                    mime="text/csv"
                )
        
        with col2:
            if 'validation' in st.session_state.processed_data:
                csv_valid = st.session_state.processed_data['validation'].to_csv(index=False)
                st.download_button(
                    label="下载验证集结果 (CSV)",
                    data=csv_valid,
                    file_name="peak_intensity_validation.csv",
                    mime="text/csv"
                )

with tab3:
    st.header("数据可视化")
    
    if st.session_state.processed_data:
        # 质谱图
        if 'spectrum' in st.session_state.processed_data:
            st.subheader("平均质谱图")
            fig_spectrum = plot_spectrum(st.session_state.processed_data['spectrum'])
            st.plotly_chart(fig_spectrum, use_container_width=True)
        
        # 热图
        if 'train' in st.session_state.processed_data:
            st.subheader("训练集峰强度热图")
            fig_heatmap = plot_heatmap(st.session_state.processed_data['train'])
            st.plotly_chart(fig_heatmap, use_container_width=True)
        
        # 数据预览
        st.subheader("数据预览")
        
        data_view = st.selectbox(
            "选择数据集",
            ["训练集", "验证集"] if 'validation' in st.session_state.processed_data else ["训练集"]
        )
        
        if data_view == "训练集":
            st.dataframe(st.session_state.processed_data['train'], use_container_width=True)
        else:
            st.dataframe(st.session_state.processed_data['validation'], use_container_width=True)
    
    else:
        st.info("💡 请先在「处理与结果」页面处理数据")

# 页脚
st.divider()
st.markdown("""
<div style='text-align: center; color: #888; padding: 2rem 0;'>
    <p>MALDI-TOF MS 数据处理平台 | Powered by Streamlit & MALDIquant</p>
</div>
""", unsafe_allow_html=True)
