from structman.lib.output import out_generator

def init_feature_table(feature_file, obj_only = False):
    feature_output = out_generator.OutputGenerator()

    if obj_only:
        return feature_output

    feat_f = open(feature_file, 'a')
    return feature_output, feat_f
