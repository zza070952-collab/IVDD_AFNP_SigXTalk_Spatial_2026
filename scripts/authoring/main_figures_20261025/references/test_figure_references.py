import unittest
from figure_references import rewrite, scientific_numbers


class ExactMainFigureMapping(unittest.TestCase):
    def test_single_data_panels(self):
        for old, new in [('1B','2A'),('1E','2D'),('2D','3D'),('3F','4F'),('4A','5A'),('5B','6A'),('5G','6F')]:
            self.assertEqual(rewrite('Figure '+old)[0], 'Figure '+new)

    def test_deleted_workflow_semantics(self):
        self.assertEqual(rewrite('Figure 1A')[0], 'Figure 1')
        self.assertEqual(rewrite('Figure 5A')[0], 'Figure 1E')

    def test_range_and_multiple_figures(self):
        self.assertEqual(rewrite('Figures 1B–E and 5B–G')[0], 'Figures 2A–D and 6A–F')
        self.assertEqual(rewrite('Fig. 2C–F')[0], 'Figure 3C–F')
        self.assertEqual(rewrite('Figures 2 and 3')[0], 'Figures 3 and 4')

    def test_supplementary_excluded(self):
        text='Figure S1A and Supplementary Figures S2–S4 remain 0.034, n=8.'
        self.assertEqual(rewrite(text)[0], text)

    def test_no_numeric_or_citation_change(self):
        before='8 donors [1–3] (Figure 1B–E); r=0.82, P=1.4e-5, 48 μm; Figure S2.'
        after,_=rewrite(before)
        self.assertEqual(scientific_numbers(before), scientific_numbers(after))

    def test_no_recursive_replacement(self):
        self.assertEqual(rewrite('Figure 2; Figure 3; Figure 4; Figure 5')[0],
                         'Figure 3; Figure 4; Figure 5; Figure 6')

    def test_caption_separator_spacing(self):
        self.assertEqual(rewrite('Figure 1 | Caption.')[0],'Figure 2 | Caption.')


if __name__=='__main__':
    unittest.main()
