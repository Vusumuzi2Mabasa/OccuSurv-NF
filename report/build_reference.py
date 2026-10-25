"""Rebuild the committed Quarto Word style reference with python-docx."""
from pathlib import Path
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.style import WD_STYLE_TYPE
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
d=Document();s=d.sections[0]
s.page_width=Inches(8.5);s.page_height=Inches(11)
s.top_margin=s.bottom_margin=Inches(.7);s.left_margin=s.right_margin=Inches(.8)
for st in d.styles:
    for border in list(st.element.iter(qn('w:pBdr'))):border.getparent().remove(border)
for name in ['Normal','Body Text','Title','Subtitle','Heading 1','Heading 2','Heading 3','Caption']:
    st=d.styles[name];st.font.name='Calibri';st.font.color.rgb=RGBColor(0,0,0);st.font.size=Pt(11)
    st.paragraph_format.space_after=Pt(7)
for name,size in [('Title',25),('Subtitle',14),('Heading 1',17),('Heading 2',15),('Heading 3',12),('Caption',9)]:d.styles[name].font.size=Pt(size)
st=d.styles.add_style('Table',WD_STYLE_TYPE.TABLE);st.font.name='Calibri';st.font.size=Pt(10)
pr=OxmlElement('w:tblPr');b=OxmlElement('w:tblBorders')
for edge in ['top','left','bottom','right','insideH','insideV']:
    x=OxmlElement('w:'+edge)
    for key,val in [('val','single'),('sz','4'),('color','D9D9D9')]:x.set(qn('w:'+key),val)
    b.append(x)
pr.append(b);m=OxmlElement('w:tblCellMar')
for edge in ['top','left','bottom','right']:
    x=OxmlElement('w:'+edge);x.set(qn('w:w'),'60');x.set(qn('w:type'),'dxa');m.append(x)
pr.append(m);st.element.append(pr)
ppr=OxmlElement('w:pPr');sp=OxmlElement('w:spacing');sp.set(qn('w:after'),'40');ppr.append(sp);st.element.append(ppr)
first=OxmlElement('w:tblStylePr');first.set(qn('w:type'),'firstRow')
rp=OxmlElement('w:rPr');rp.append(OxmlElement('w:b'));first.append(rp)
tc=OxmlElement('w:tcPr');shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'E7EEF3');tc.append(shade);first.append(tc);st.element.append(first)
d.add_paragraph('OccuSurv reference styles',style='Title')
d.save(Path(__file__).with_name('reference.docx'))
