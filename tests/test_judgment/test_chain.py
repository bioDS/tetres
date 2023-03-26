import pytest
import random
from pathlib import Path
from tetres.judgment.chain import Chain
import os


def test_MChain_construction_files():
    assert Chain(
        trees="30Taxa.trees",
        log_file="30Taxa_beast2.log",
        summary="30Taxa_centroid.tree",
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
    ), "Construction files failed"


def test_MChain_construction_TimeTreeSets(thirty_taxa_tts, thirty_taxa_centroid_tts):
    assert Chain(
        trees=thirty_taxa_tts,
        log_file="30Taxa_beast2.log",
        summary=thirty_taxa_centroid_tts,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
    ), "Construction TTS failed"


def test_MChain_construction_mixed1(thirty_taxa_tts):
    assert Chain(
        trees=thirty_taxa_tts,
        log_file="30Taxa_beast2.log",
        summary="30Taxa_centroid.tree",
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
    ), "Construction TTS trees, centroid file failed"


def test_MChain_construction_mixed2(thirty_taxa_centroid_tts):
    assert Chain(
        trees="30Taxa.trees",
        log_file="30Taxa_beast2.log",
        summary=thirty_taxa_centroid_tts,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
    ), "Construction trees file, TTS centroid failed"


def test_MChain_construction_summary_None():
    assert Chain(
        trees="30Taxa.trees",
        log_file="30Taxa_beast2.log",
        summary=None,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
    ), "Construction with summary=None failed"


def test_MChain_construction_workingdir_None():
    with pytest.raises(ValueError):
        Chain(
            trees="30Taxa.trees",
            log_file="30Taxa_beast2.log",
            summary="30Taxa_centroid.tree",
            working_dir=None
        )


# testing construction exceptions
def test_MChain_wrong_trees():
    with pytest.raises(ValueError):
        Chain(
            trees=20,
            log_file="30Taxa_beast2.log",
            summary="30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
        )


def test_MChain_wrong_logfile_path():
    with pytest.raises(FileNotFoundError):
        Chain(
            trees="30Taxa.trees",
            log_file="30Taxa.log",
            summary="30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
        )


def test_MChain_wrong_logfile_type():
    with pytest.raises(ValueError):
        Chain(
            trees="30Taxa.trees",
            log_file=20,
            summary="30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
        )


def test_MChain_wrong_summary():
    with pytest.raises(ValueError):
        Chain(
            trees="30Taxa.trees",
            log_file="30Taxa_beast2.log",
            summary=20,
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
        )


def test_MChain_wrong_workingdir_path():
    with pytest.raises(NotADirectoryError):
        Chain(
            trees="30Taxa.trees",
            log_file="30Taxa_beast2.log",
            summary="30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/work"
        )


def test_MChain_wrong_workingdir_type():
    with pytest.raises(ValueError):
        Chain(
            trees="30Taxa.trees",
            log_file="30Taxa_beast2.log",
            summary="30Taxa_centroid.tree",
            working_dir=20
        )


def test_MChain_wrong_treemaps():
    with pytest.raises(ValueError):
        Chain(
            trees="30Taxa.trees",
            log_file="30Taxa_beast2.log",
            summary="12Taxa.trees",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
        )


def test_MChain_get_key_names(thirty_taxa_chain):
    assert thirty_taxa_chain.get_key_names() == ['posterior',
                                                  'likelihood',
                                                  'prior',
                                                  'treeLikelihood',
                                                  'TreeHeight',
                                                  'popSize',
                                                  'CoalescentConstant'], \
        "Get_key_name functions failed!"


def test_MChain_get_ess_partial_tracerer(thirty_taxa_chain):
    result = []
    for ess_key in thirty_taxa_chain.get_key_names():
        result.append(thirty_taxa_chain.get_ess(ess_key=ess_key, upper_i=200, lower_i=100))
    assert result == [99.01960784313721,
                         99.01960784313725,
                         99.01960784313725,
                         99.01960784313725,
                         89.81794714804302,
                         99.01960784313725,
                         99.01960784313725], \
        "ESS partial value calculation with MChain failed!"


def test_MChain_get_ess_partial_wronglowertype(thirty_taxa_chain):
    with pytest.raises(ValueError):
        thirty_taxa_chain.get_ess(ess_key="posterior", lower_i=10.2)


def test_MChain_get_ess_partial_wronguppertype(thirty_taxa_chain):
    with pytest.raises(ValueError):
        thirty_taxa_chain.get_ess(ess_key="posterior", upper_i=10.2)


def test_MChain_get_ess_partial_uppertoobig(thirty_taxa_chain):
    with pytest.raises(IndexError):
        thirty_taxa_chain.get_ess(ess_key="posterior", upper_i=3000)


def test_MChain_get_ess_partial_wrongway(thirty_taxa_chain):
    with pytest.raises(ValueError):
        thirty_taxa_chain.get_ess(ess_key="posterior", lower_i=10, upper_i=5)


def test_MChain_get_ess_partial_neglower(thirty_taxa_chain):
    with pytest.raises(ValueError):
        thirty_taxa_chain.get_ess(ess_key="posterior", lower_i=-10)


def test_MChain_get_ess_partial_negupper(thirty_taxa_chain):
    with pytest.raises(ValueError):
        thirty_taxa_chain.get_ess(ess_key="posterior", upper_i=-10)


def test_MChain_get_ess_partial_correctness(thirty_taxa_chain):
    result = []
    for ess_key in thirty_taxa_chain.get_key_names():
        result.append(thirty_taxa_chain.get_ess(ess_key=ess_key))

    result_partial = []
    for ess_key in thirty_taxa_chain.get_key_names():
        result_partial.append(thirty_taxa_chain.get_ess(ess_key=ess_key, lower_i=0,
                                                        upper_i=thirty_taxa_chain.log_data.shape[0] - 1))
    assert result == result_partial, \
        "ESS partial returns wrong ESS!"


def test_MChain_pseudo_ess(thirty_taxa_chain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    p_ess = thirty_taxa_chain.get_pseudo_ess(sample_range=len(thirty_taxa_chain))
    random.setstate(state)  # reset the random seed to previous state
    assert int(p_ess) == 897, "Get Pseudo ESS for MChain failed!"


def test_MChain_pseudo_ess_few_samples(thirty_taxa_chain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    p_ess = thirty_taxa_chain.get_pseudo_ess(sample_range=25)
    random.setstate(state)  # reset the random seed to previous state
    assert int(p_ess) == 896, "Get Pseudo ESS for MChain failed!"


def test_MChain_pwd_matrix(thirty_taxa_chain):
    thirty_taxa_chain.pwd_matrix()
    assert os.path.exists(f"{thirty_taxa_chain.working_dir}/data/30TestFix.npy"),\
            "Similarity matrix computation failed!"


def test_MChain_pwd_matrix(thirty_taxa_chain):
    pwd = thirty_taxa_chain.pwd_matrix(rf=True)
    assert os.path.exists(f"{thirty_taxa_chain.working_dir}/data/30TestFix_rf.npy"),\
            "Similarity matrix computation failed!"


def test_MChain_get_simmatrix(thirty_taxa_chain):
    for beta in [1, 2, 0.5]:
        thirty_taxa_chain.simmilarity_matrix(beta=beta)
        assert os.path.exists(f"{thirty_taxa_chain.working_dir}/data/30TestFix_{'' if beta == 1 else f'{beta}_'}similarity.npy"),\
            "Similarity matrix computation failed!"


# Currently not supported feature, WIP
# def test_MChain_spectral_clustree(thirty_taxa_chain):
#     for beta in [1, 2, 0.5]:
#         thirty_taxa_chain.spectral_clustree(beta=beta)
#         assert os.path.exists(f"{thirty_taxa_chain.working_dir}/data/30TestFix_{'' if beta == 1 else f'{beta}_'}clustering.npy"),\
#             "Spectral clustering for mchain failed!"


def test_MChain_split_trees_from_clustering(ten_taxa_multichain):
    # todo throws error for already existing tree files
    ten_taxa_multichain[0].split_trees_from_clustering(k=2)
