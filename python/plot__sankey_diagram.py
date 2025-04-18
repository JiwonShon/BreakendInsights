import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = 'browser'


df = similarity_df.copy()
patientID='patientid'

df['patient_barcode'] = df['file_name'].str.extract(r'^(.*)_amplicon')
df['aliquot_barcode1'] = df['Amp1'].str.extract(r'^(.*)_amplicon')
df['aliquot_barcode2'] = df['Amp2'].str.extract(r'^(.*)_amplicon')
df['tumor_stage1'] = df['aliquot_barcode1'].str.split('__').str[0].str[12:]
df['tumor_stage2'] = df['aliquot_barcode2'].str.split('__').str[0].str[12:]

# 그다음 amplicon, ampType, amp_num 순으로 추출
df[['amplicon1']] = df['Amp1'].str.extract(r'_(amplicon\d+)')
df[['amplicon2']] = df['Amp2'].str.extract(r'_(amplicon\d+)')

df = df.merge(categories_df.rename(columns={'Amp': 'Amp1', 'CAT': 'CAT1'}), on='Amp1', how='left')
df = df.merge(categories_df.rename(columns={'Amp': 'Amp2', 'CAT': 'CAT2'}), on='Amp2', how='left')

test=df[df['patient_barcode']==patientID][['patient_barcode','tumor_stage1','tumor_stage2','amplicon1','amplicon2','CAT1','CAT2','SimilarityScore','SimScorePercentile']]

df = test.copy()
print(test.shape)

# -------------------- 설정 --------------------
stage_order = ['T', 'TII', 'TIII', 'TIV']
stage_rank = {stage: i for i, stage in enumerate(stage_order)}

cat_color_map = {
    'ecDNA+': 'red',
    'ecDNA-': 'orange',
    'Linear amplification': 'blue',
    'Complex non-cyclic': 'yellow',
    'No amp/Invalid': 'gray',
    None: 'lightgray'
}

# -------------------- 데이터 정리 --------------------

# Stage 순서 정렬 (낮은 쪽을 stage1으로 고정)
def sort_stages(row):
    r1 = stage_rank.get(row['tumor_stage1'], -1)
    r2 = stage_rank.get(row['tumor_stage2'], -1)
    if r1 > r2:
        return pd.Series({
            "tumor_stage1": row['tumor_stage2'],
            "tumor_stage2": row['tumor_stage1'],
            "amplicon1": row['amplicon2'],
            "amplicon2": row['amplicon1'],
            "CAT1": row['CAT2'],
            "CAT2": row['CAT1'],
            "SimilarityScore": row['SimilarityScore']
        })
    else:
        return row[["tumor_stage1", "tumor_stage2", "amplicon1", "amplicon2", "CAT1", "CAT2", "SimilarityScore"]]

df = df.apply(sort_stages, axis=1)

# stage rank 매기고 연속적 transition만 필터
df['stage1_rank'] = df['tumor_stage1'].map(stage_rank)
df['stage2_rank'] = df['tumor_stage2'].map(stage_rank)
df = df[df['stage2_rank'] - df['stage1_rank'] == 1].copy()

# source/target 노드와 값 지정
df['source'] = df['tumor_stage1'] + '_' + df['amplicon1']
df['target'] = df['tumor_stage2'] + '_' + df['amplicon2']
df['value'] = df['SimilarityScore']
df['source_cat'] = df['CAT1']

# -------------------- 노드 정보 정리 --------------------
node_cat_df = pd.concat([
    df[['source', 'source_cat']].rename(columns={'source': 'node', 'source_cat': 'cat'}),
    df[['target', 'source_cat']].rename(columns={'target': 'node', 'source_cat': 'cat'})
]).drop_duplicates('node')

def get_stage(node_name):
    return node_name.split('_')[0]

# 노드 정렬 및 색상 매핑
all_nodes_sorted = sorted(node_cat_df['node'].unique(), key=lambda x: (stage_rank.get(get_stage(x), 999), x))
node_cat_df = node_cat_df.set_index('node').reindex(all_nodes_sorted)
node_colors = node_cat_df['cat'].map(cat_color_map).fillna('lightgray').tolist()

# 노드 인덱싱
node_dict = {name: idx for idx, name in enumerate(all_nodes_sorted)}
df['source_id'] = df['source'].map(node_dict)
df['target_id'] = df['target'].map(node_dict)

# -------------------- 좌표 배치 --------------------
node_x = [stage_rank.get(get_stage(n), 0) / (len(stage_order)-1) for n in all_nodes_sorted]

# 세로 위치 압축 배치
spacing = 0.02
node_y = np.arange(0, len(all_nodes_sorted)) * spacing
if node_y.max() > 1.0:
    node_y = node_y / node_y.max()

# -------------------- [추가] 흐름선 색상 진하기 조절 --------------------
score_min = df['value'].min()
score_max = df['value'].max()
df['alpha'] = df['value'].apply(lambda v: 0.2 + 0.8 * (v - score_min) / (score_max - score_min))
df['link_color'] = df['alpha'].apply(lambda a: f'rgba(80,80,80,{a:.2f})')
print(df.shape)

# -------------------- Sankey 그리기 --------------------
fig = go.Figure(data=[go.Sankey(
    arrangement='snap',
    node=dict(
        pad=5,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=all_nodes_sorted,
        color=node_colors,
        x=node_x,
        y=node_y
    ),
    link=dict(
        source=df['source_id'],
        target=df['target_id'],
        value=df['value'],
        color=df['link_color']  # ✅ 진하기 반영된 회색 링크
    )
)])

# -------------------- 링크 hover 텍스트 설정 --------------------
link_labels = df.apply(
    lambda row: f"{row['source']} → {row['target']}<br>SimilarityScore: {row['value']:.3f}",
    axis=1
)
fig['data'][0]['link']['customdata'] = link_labels
fig['data'][0]['link']['hovertemplate'] = '%{customdata}<extra></extra>'

# -------------------- 범례용 dummy trace 삽입 --------------------
for cat, color in cat_color_map.items():
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color=color),
        name=cat,
        legendgroup=cat,
        showlegend=True
    ))

# -------------------- 전체 레이아웃 업데이트 --------------------
fig.update_layout(
    title_text=f"Sankey Diagram: Adjacent Tumor Stage Transitions Only: {patientID}",
    font_size=10,
    width=1200,
    height=600,
    paper_bgcolor='rgba(0,0,0,0)',   # ✅ 배경 제거
    plot_bgcolor='rgba(0,0,0,0)',    # ✅ 배경 제거
    legend=dict(
        title="CAT (Node Color)",
        orientation="v",
        x=1.02,
        y=1,
        xanchor="left",
        yanchor="top"
    ),
    xaxis=dict(visible=False),  # ✅ x축 제거
    yaxis=dict(visible=False)   # ✅ y축 제거
)

fig.show()
