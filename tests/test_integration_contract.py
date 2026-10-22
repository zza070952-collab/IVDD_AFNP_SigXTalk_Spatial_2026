"""Synthetic structure-only tests; no fixture is research evidence."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/09_integrate_evidence.py'
spec=importlib.util.spec_from_file_location('integration',SCRIPT)
integration=importlib.util.module_from_spec(spec)
spec.loader.exec_module(integration)

def technical_b():
    full=pd.DataFrame([
        dict(Receptor='R',TF='SSC',TG='PLOD2',pred_label=.9,PRS=.2,PRS_relative_max_filter_pass=True),
        dict(Receptor='R',TF='SSC',TG='CCN2',pred_label=.2,PRS=np.nan,PRS_relative_max_filter_pass=False),
        dict(Receptor='R',TF='SSC',TG='THBS1',pred_label=.9,PRS=np.nan,PRS_relative_max_filter_pass=False),
        dict(Receptor='R',TF='SSC',TG='PRSS23',pred_label=.9,PRS=.001,PRS_relative_max_filter_pass=False)])
    context=pd.DataFrame([dict(source='AF_stromal',target='NP_GEP2_high',ligand='FN1',
        receptor='ITGA3_ITGB1',interaction_name='TECHNICAL_FIXTURE_ONLY',Receptor='R',TF='SSC',TG='PLOD2',
        pred_label=.9,PRS=.2,prob=.01,axis_key='R::SSC::PLOD2',axis_source_row=2,
        lr_source_row=2,complex_contains_ITGB1=True)])
    return dict(full=full,context=context,activation_threshold=.75,testability=pd.DataFrame(),
                selection=pd.DataFrame({'gene':['PLOD2','CCN2','THBS1','PRSS23','LOX'],
                                        'selected_target':[True,True,True,True,False]}))

class IntegrationContractTest(unittest.TestCase):
    def test_partial_document_generation_is_explicitly_pending(self):
        # Entire fixture lives in a temporary directory and is deleted. Its
        # numeric values test file plumbing only, never biological conclusions.
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);a=root/'01_human_programs';c=root/'03_spatial_gse287499'
            a.mkdir();c.mkdir()
            (a/'STAGE_A_COMPUTATION_STATUS.json').write_text(json.dumps(dict(
                Tu_donors=8,Swahn_donors=8,Tu_cells_qc_pass=80,Tu_score_method='TECHNICAL_FIXTURE',
                Swahn_score_origin='TECHNICAL_FIXTURE_NOT_RESEARCH')))
            pd.DataFrame([dict(program='GEP2',rank=i+1,gene=g) for i,g in enumerate(integration.CANDIDATES)]).to_csv(
                a/'FROZEN_K5_PROGRAMS.tsv',sep='\t',index=False)
            pd.DataFrame([dict(program='GEP2',mapped_genes_min=10,mapped_genes_max=10,
                NP_donor_median_min=.1,NP_donor_median_max=.2,expression_interpretation='TECHNICAL_FIXTURE_NOT_RESEARCH')]).to_csv(
                a/'ALL_PROGRAMS_BRIEF.tsv',sep='\t',index=False)
            pd.DataFrame([dict(cohort='Tu_TECHNICAL',donor_id=str(i),
                NP_gt_context=True if i<7 else np.nan,
                context_donor_median=.1 if i<7 else np.nan) for i in range(8)]).to_csv(
                a/'GEP2_DONOR_SUPPORT.tsv',sep='\t',index=False)
            pd.DataFrame({'fixture':['NOT_RESEARCH']}).to_csv(a/'PROGRAM_SCORE_DONOR_SUMMARY.tsv',sep='\t',index=False)
            tu=[];mouse=[]
            for gene in integration.CANDIDATES:
                for i in range(8):
                    for identity,val in [('NP_matrix_marker_supported',.2),('immune_marker_supported',.1)]:
                        tu.append(dict(donor_id='TECHNICAL_'+str(i),identity=identity,gene=gene,n_cells=10,
                            fraction_detected=val,mean_log1p_published_expression=val,feature_present=True,
                            evidence_type='SYNTHETIC_TECHNICAL_FIXTURE_NOT_RESEARCH'))
                    for region in ['NP_CORE','AF_CORE']:
                        mouse.append(dict(human_gene=gene,capture_id='FIXTURE_'+str(i//4),disc_id='disc'+str(i),
                            region=region,n_bins=10,fraction_detected_bins=.1,n_detected_bins=1,
                            total_raw_UMI=1,native_feature_assessable=True))
            pd.DataFrame(tu).to_csv(a/'TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv',sep='\t',index=False)
            pd.DataFrame(mouse).to_csv(c/'CANDIDATE_GENES_PER_DISC_REGION.tsv',sep='\t',index=False)
            (c/'DEPTH_ANALYSIS_STATUS.json').write_text(json.dumps(dict(
                status='DONE_SINGLE_PREDECLARED_DEPTH_CHECK',protected_old_inputs_unchanged=True,
                outcome='TECHNICAL_FIXTURE_NOT_RESEARCH',adequate_paired_discs=7,
                NP_greater_AF_after_depth100=6,mean_paired_delta_original_same_eligible_bins=.2,
                mean_paired_delta_depth100=.1,observed_mean_difference_attenuated=True)))
            result=subprocess.run([sys.executable,str(SCRIPT),'--project-root',str(root),'--allow-incomplete-b'],
                                  capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            status=json.loads((root/'04_integrated/INTEGRATION_STATUS.json').read_text())
            self.assertEqual(status['status'],'PARTIAL_NEW_B_PENDING')
            self.assertFalse(status['formal_new_B_verified'])
            self.assertEqual(status['github_code_release'],'BLOCKED_UNVERIFIED')
            results=(root/'06_manuscript/RESULTS_COMPUTATION_ACTUALLY_EXECUTED.md').read_text(encoding='utf8')
            self.assertIn('No new formal pathway count',results)
            self.assertIn('7/7 assessable Tu donors (of 8 total)',results)
            self.assertFalse((root/'04_integrated/CANDIDATE_ALL_TESTED_NEW_AXES_WITH_PRS.tsv.gz').exists())
            # Continue the disposable structure fixture through the completion
            # gate to exercise joins/report generation. These are simulated
            # receipts in a deleted temp directory, not execution evidence.
            wd=root/'02_sigxtalk_official_new/formal_new'
            for folder in ['audit','tables','prs_outputs']:(wd/folder).mkdir(parents=True,exist_ok=True)
            prov=dict(input_kind='real_GSE230808_raw_counts',started_from_old_pathway_output=False,
                      project_id='TECHNICAL_FIXTURE_NOT_RESEARCH',prepare_parameters={'CCC_threshold':.1})
            (wd/'audit/NEW_INPUT_PROVENANCE.json').write_text(json.dumps(prov))
            comp=dict(complete_R_python_R=True,no_old_pathways_or_PRS_used=True,
                new_HGNN_predictions_sha256='TECHNICAL_HASH',
                new_input_provenance_sha256=integration.sha(wd/'audit/NEW_INPUT_PROVENANCE.json'),
                all_tested_HGNN_rows=4,positive_PRS_rows=2,active_PA_rows=3,filtered_PRS_rows=1,
                fidelity_specificity_status='DONE')
            (wd/'audit/NEW_R_PYTHON_R_COMPLETION.json').write_text(json.dumps(comp))
            (wd/'NEW_OFFICIAL_FINAL_STAGE_STATUS.json').write_text(json.dumps(dict(
                formal_design=True,new_result_not_legacy_relabel=True,run_class='FULL_NEW')))
            (wd/'audit/HGNN_STATUS.json').write_text(json.dumps(dict(status='DONE',new_computation=True,
                mode='formal',prediction_sha256='TECHNICAL_HASH',epochs_completed=100,seed=123)))
            full=technical_b()['full'].copy();full['Receptor']='ITGB1'
            full['method_origin']='NEW_OFFICIAL_SIGXTALK_TECHNICAL_FIXTURE'
            full['source_project_id']=prov['project_id']
            full['source_run_started_from_legacy_predictions']=False
            full.to_csv(wd/'tables/ALL_TESTED_NEW_AXES_WITH_PRS.tsv',sep='\t',index=False)
            lr=dict(source='AF_stromal',target='NP_GEP2_high',ligand='FN1',receptor='ITGA3_ITGB1',
                    interaction_name='TECHNICAL_FIXTURE_ONLY',prob=.01,pval=.01)
            pd.DataFrame([lr]).to_csv(wd/'tables/FOCUS_LR_COMPLEX_ALL.tsv',sep='\t',index=False)
            components=[]
            for component in ['ITGA3','ITGB1']:
                z={k:v for k,v in lr.items() if k not in ['prob','pval','receptor']}
                z.update(receptor_complex='ITGA3_ITGB1',receptor_component=component);components.append(z)
            pd.DataFrame(components).to_csv(wd/'tables/LR_COMPLEX_COMPONENTS.tsv',sep='\t',index=False)
            pd.DataFrame([dict(From='FN1',To='ITGB1',Weight=.01,ligand_receiver_gene_present=True,
                              receptor_receiver_gene_present=True)]).to_csv(
                wd/'tables/PREPARE_LR_RECEIVER_FILTER_AUDIT.tsv',sep='\t',index=False)
            pd.DataFrame([dict(Receptor='ITGB1',SSC='SSC',Target='PLOD2',PRS=.2,
                fidelity_target_all=.2,specificity_receptor_all=.3,normalization_universe='TECHNICAL_ONLY')]).to_csv(
                wd/'prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv',sep='\t',index=False)
            formal=subprocess.run([sys.executable,str(SCRIPT),'--project-root',str(root)],capture_output=True,text=True)
            self.assertEqual(formal.returncode,0,formal.stdout+formal.stderr)
            formal_status=json.loads((root/'04_integrated/INTEGRATION_STATUS.json').read_text())
            self.assertEqual(formal_status['status'],'DONE_NEW_FORMAL_INTEGRATION')
            self.assertTrue(formal_status['FN1_ITGB1_complex_context_with_PLOD2'])
            matrix=pd.read_csv(root/'04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv',sep='\t').set_index('gene')
            self.assertEqual(matrix.loc['PLOD2','n_axes_with_official_filtered_universe_fid_spe'],1)

    def test_missing_formal_status_is_rejected_without_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            result=subprocess.run([sys.executable,str(SCRIPT),'--project-root',td],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('New formal SigXTalk completion is not verified',result.stderr)
            self.assertFalse((Path(td)/'04_integrated').exists())

    def test_target_reasons_keep_NA_and_threshold_failures_separate(self):
        b=technical_b()
        expected={'PLOD2':'SOURCE_SUPPORTED_RECEPTOR_KEY_CANDIDATE_CONTEXT',
                  'CCN2':'TESTED_NO_PA_THRESHOLD_PASS',
                  'THBS1':'PA_PASS_NO_RETURNED_OFFICIAL_PRS',
                  'PRSS23':'POSITIVE_PRS_NO_RELATIVE_FILTER_PASS',
                  'LOX':'NOT_IN_TESTED_TARGET_UNIVERSE'}
        for gene,status in expected.items():
            result=integration.target_evidence(gene,b,True)
            self.assertEqual(result['target_testability'],status)
        self.assertTrue(np.isnan(integration.target_evidence('THBS1',b,True)['max_original_PRS']))
        details=json.loads(integration.target_evidence('PLOD2',b,True)['representative_path_details_json'])
        self.assertEqual(details[0]['axis_key'],'R::SSC::PLOD2')
        self.assertEqual(details[0]['PRS'],.2)

    def test_unfinished_B_cannot_be_rescued_by_available_rows(self):
        result=integration.target_evidence('PLOD2',technical_b(),False)
        self.assertEqual(result['target_testability'],'NOT_ASSESSED_NEW_FORMAL_RUN_INCOMPLETE')
        self.assertTrue(np.isnan(result['n_tested_axes']))
        self.assertEqual(result['representative_path_details_json'],'[]')

    def test_code_release_requires_explicit_verified_public_receipt(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'proof.json'
            path.write_text(json.dumps({'repo_url':'https://github.com/example/repo','head_sha':'a'*40}))
            self.assertFalse(integration.verified_release(path)[0])
            path.write_text(json.dumps({'verified_public':True,'repo_url':'https://github.com/example/repo',
                                        'head_sha':'abc','verified_utc':'TEST_FIXTURE'}))
            self.assertFalse(integration.verified_release(path)[0])

if __name__=='__main__': unittest.main()
