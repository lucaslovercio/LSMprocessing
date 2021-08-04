from keras import backend as K
import tensorflow as tf
import numpy as np
# K.epsilon returns 1e-7

#TODO
'''
- debug why the epochwise metrics aren't working on compute canada
- debug why there's no difference b/w training and validation F1 when computed on whole epochs
 - http://digital-thinking.de/keras-three-ways-to-use-custom-validation-metrics-in-keras/
 - https://keras.io/guides/writing_your_own_callbacks/
- F1macro etc relies on there being exactly 3 classes, should be made to generalize better probably

'''

'''
custom metrics:
need to take y_true and y_pred as arguments
need to return a single tensor value
the returned tensors have shape (# of rows, 1)

y_pred:
The output of the last model layer, softmax probabilities in a vector of size <# of classes> for the current model.
Has shape = (batch_size, 128, 128, 3).

y_true:
The ground truth segmentations, 1-hot encoded.
Has shape = (batch_size, 128, 128, 3).
'''

########### below here are stateless versions of metrics

# compute the confusion matrix values for a given class
def get_confusion_matrix_values(y_true, y_pred, class_num):
    # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
    y_pred_flat = tf.math.argmax(y_pred, axis = -1)
    y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
    # slice y_true and y_pred to only consider one class
    y_pred = y_pred[:,:,:,class_num]
    y_true = y_true[:,:,:,class_num]
    # construct matrix represention of where predictions matched the ground truth
    correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
    # construct matrix represention of where predictions didn't match the ground truth
    incorrect = 1.0 - correct
    # compute the true and false positives and negatives
    tp = tf.reduce_sum(y_pred * correct)
    tn = tf.reduce_sum((1.0 - y_pred) * correct)
    fp = tf.reduce_sum(y_pred * incorrect)
    fn = tf.reduce_sum((1.0 - y_pred) * incorrect)
    return tp, tn, fp, fn

# compute the recall for a given class
def recall_class(y_true, y_pred, class_num):
    tp, _, _, fn = get_confusion_matrix_values(y_true, y_pred, class_num)
    return tp / (tp + fn + K.epsilon())

# compute the precision for a given class
def precision_class(y_true, y_pred, class_num):
    tp, _, fp, _ = get_confusion_matrix_values(y_true, y_pred, class_num)
    return tp / (tp + fp + K.epsilon())

# macroaveraged recall over all classes
def recall_macro_batch(y_true, y_pred):
    num_classes = y_pred.shape[-1]
    recall_accumulator = 0
    for i in range(num_classes):
        recall_accumulator += recall_class(y_true, y_pred, i)
    return recall_accumulator / num_classes

# macroaveraged precision over all classes
def precision_macro_batch(y_true, y_pred):
    num_classes = y_pred.shape[-1]
    precision_accumulator = 0
    for i in range(num_classes):
        precision_accumulator += precision_class(y_true, y_pred, i)
    return precision_accumulator / num_classes

# macroaveraged f1 score over all classes
def f1_macro_batch(y_true, y_pred):
    recall = recall_macro_batch(y_true, y_pred)
    precision = precision_macro_batch(y_true, y_pred)
    return 2 * (precision * recall) / (precision + recall + K.epsilon())

# precision and recall for individual classes
def prec_0(y_true, y_pred):
    return precision_class(y_true, y_pred, 0)

def reca_0(y_true, y_pred):
    return recall_class(y_true, y_pred, 0)

def prec_1(y_true, y_pred):
    return precision_class(y_true, y_pred, 1)

def reca_1(y_true, y_pred):
    return recall_class(y_true, y_pred, 1)

def prec_2(y_true, y_pred):
    return precision_class(y_true, y_pred, 2)

def reca_2(y_true, y_pred):
    return recall_class(y_true, y_pred, 2)


# accuracy b/w the two types of tissues
def tissue_type_accuracy_batch(y_true, y_pred):
    y_true = tf.argmax(y_true, axis=-1)
    y_pred = tf.argmax(y_pred, axis=-1)

    tissue_indices = tf.cast(y_true, dtype=tf.bool)

    y_pred = tf.boolean_mask(y_pred, tissue_indices)
    y_true = tf.boolean_mask(y_true, tissue_indices)

    total_correct = tf.reduce_sum(tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32))
    
    return total_correct / tf.cast(tf.size(y_pred), dtype=tf.float32)


# accuracy of tissue vs. background segmentations
def binary_accuracy_batch(y_true, y_pred):
    y_pred_flat = tf.math.argmax(y_pred, axis = -1)
    y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
    # remove the predictions and ground truth for the non-background classes
    y_pred = y_pred[:,:,:,0]
    y_true = y_true[:,:,:,0]
    # construct matrix represention of where predictions matched the ground truth
    correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
    return tf.reduce_sum(correct) / tf.cast(tf.size(correct), dtype=tf.float32)

########### below here are stateful versions of metrics

# tensorflow documentation https://www.tensorflow.org/api_docs/python/tf/keras/metrics/Metric
# other documentation https://neptune.ai/blog/keras-metrics

class F1Macro(tf.keras.metrics.Metric):
    def __init__(self, name='f1_macro', **kwargs):
        super(F1Macro, self).__init__(name=name, **kwargs)
        self.true_positives_0 = self.add_weight(name='tp0', initializer='zeros')
        self.false_positives_0 = self.add_weight(name='fp0', initializer='zeros')
        self.false_negatives_0 = self.add_weight(name='fn0', initializer='zeros')

        self.true_positives_1 = self.add_weight(name='tp1', initializer='zeros')
        self.false_positives_1 = self.add_weight(name='fp1', initializer='zeros')
        self.false_negatives_1 = self.add_weight(name='fn1', initializer='zeros')

        self.true_positives_2 = self.add_weight(name='tp2', initializer='zeros')
        self.false_positives_2 = self.add_weight(name='fp2', initializer='zeros')
        self.false_negatives_2 = self.add_weight(name='fn2', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives_0.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * correct[:,:,:,0]))
        self.false_positives_0.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * incorrect[:,:,:,0]))
        self.false_negatives_0.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,0]) * incorrect[:,:,:,0]))
        # for neural
        self.true_positives_1.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * correct[:,:,:,1]))
        self.false_positives_1.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * incorrect[:,:,:,1]))
        self.false_negatives_1.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,1]) * incorrect[:,:,:,1]))
        # for mesen
        self.true_positives_2.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * correct[:,:,:,2]))
        self.false_positives_2.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * incorrect[:,:,:,2]))
        self.false_negatives_2.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,2]) * incorrect[:,:,:,2]))

    def result(self):
        precision_0 = (K.epsilon() + self.true_positives_0) / (self.true_positives_0 + self.false_positives_0 + K.epsilon())#new
        recall_0 = (K.epsilon() + self.true_positives_0) / (self.true_positives_0 + self.false_negatives_0 + K.epsilon())#new

        precision_1 = (K.epsilon() + self.true_positives_1) / (self.true_positives_1 + self.false_positives_1 + K.epsilon())#new
        recall_1 = (K.epsilon() + self.true_positives_1) / (self.true_positives_1 + self.false_negatives_1 + K.epsilon())#new

        precision_2 = (K.epsilon() + self.true_positives_2) / (self.true_positives_2 + self.false_positives_2 + K.epsilon())#new
        recall_2 = (K.epsilon() + self.true_positives_2) / (self.true_positives_2 + self.false_negatives_2 + K.epsilon())#new

        average_precision = (precision_0 + precision_1 + precision_2) / 3.0
        average_recall = (recall_0 + recall_1 + recall_2) / 3.0

        return 2 * (average_precision * average_recall + K.epsilon()) / (average_precision + average_recall + K.epsilon())#new

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives_0.assign(0.0)
        self.false_positives_0.assign(0.0)
        self.false_negatives_0.assign(0.0)

        self.true_positives_1.assign(0.0)
        self.false_positives_1.assign(0.0)
        self.false_negatives_1.assign(0.0)

        self.true_positives_2.assign(0.0)
        self.false_positives_2.assign(0.0)
        self.false_negatives_2.assign(0.0)


class PrecisionMacro(tf.keras.metrics.Metric):
    def __init__(self, name='precision_macro', **kwargs):
        super(PrecisionMacro, self).__init__(name=name, **kwargs)
        self.true_positives_0 = self.add_weight(name='tp0', initializer='zeros')
        self.false_positives_0 = self.add_weight(name='fp0', initializer='zeros')

        self.true_positives_1 = self.add_weight(name='tp1', initializer='zeros')
        self.false_positives_1 = self.add_weight(name='fp1', initializer='zeros')

        self.true_positives_2 = self.add_weight(name='tp2', initializer='zeros')
        self.false_positives_2 = self.add_weight(name='fp2', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives_0.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * correct[:,:,:,0]))
        self.false_positives_0.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * incorrect[:,:,:,0]))
        # for neural
        self.true_positives_1.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * correct[:,:,:,1]))
        self.false_positives_1.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * incorrect[:,:,:,1]))
        # for mesen
        self.true_positives_2.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * correct[:,:,:,2]))
        self.false_positives_2.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * incorrect[:,:,:,2]))

    def result(self):
        precision_0 = (K.epsilon() + self.true_positives_0) / (self.true_positives_0 + self.false_positives_0 + K.epsilon())#new
        precision_1 = (K.epsilon() + self.true_positives_1) / (self.true_positives_1 + self.false_positives_1 + K.epsilon())#new
        precision_2 = (K.epsilon() + self.true_positives_2) / (self.true_positives_2 + self.false_positives_2 + K.epsilon())#new
        average_precision = (precision_0 + precision_1 + precision_2) / 3.0
        return average_precision

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives_0.assign(0.0)
        self.false_positives_0.assign(0.0)
        self.true_positives_1.assign(0.0)
        self.false_positives_1.assign(0.0)
        self.true_positives_2.assign(0.0)
        self.false_positives_2.assign(0.0)


class RecallMacro(tf.keras.metrics.Metric):
    def __init__(self, name='recall_macro', **kwargs):
        super(RecallMacro, self).__init__(name=name, **kwargs)
        self.true_positives_0 = self.add_weight(name='tp0', initializer='zeros')
        self.false_negatives_0 = self.add_weight(name='fn0', initializer='zeros')

        self.true_positives_1 = self.add_weight(name='tp1', initializer='zeros')
        self.false_negatives_1 = self.add_weight(name='fn1', initializer='zeros')

        self.true_positives_2 = self.add_weight(name='tp2', initializer='zeros')
        self.false_negatives_2 = self.add_weight(name='fn2', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives_0.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * correct[:,:,:,0]))
        self.false_negatives_0.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,0]) * incorrect[:,:,:,0]))
        # for neural
        self.true_positives_1.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * correct[:,:,:,1]))
        self.false_negatives_1.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,1]) * incorrect[:,:,:,1]))
        # for mesen
        self.true_positives_2.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * correct[:,:,:,2]))
        self.false_negatives_2.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,2]) * incorrect[:,:,:,2]))

    def result(self):
        recall_0 = (K.epsilon() + self.true_positives_0) / (self.true_positives_0 + self.false_negatives_0 + K.epsilon())#new
        recall_1 = (K.epsilon() + self.true_positives_1) / (self.true_positives_1 + self.false_negatives_1 + K.epsilon())#new
        recall_2 = (K.epsilon() + self.true_positives_2) / (self.true_positives_2 + self.false_negatives_2 + K.epsilon())#new
        average_recall = (recall_0 + recall_1 + recall_2) / 3.0
        return average_recall

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives_0.assign(0.0)
        self.false_negatives_0.assign(0.0)
        self.true_positives_1.assign(0.0)
        self.false_negatives_1.assign(0.0)
        self.true_positives_2.assign(0.0)
        self.false_negatives_2.assign(0.0)

# epochwise binary accuracy
class BinaryAccuracy(tf.keras.metrics.Metric):
    def __init__(self, name='binary_accuracy', **kwargs):
        super(BinaryAccuracy, self).__init__(name=name, **kwargs)
        self.correct_count = self.add_weight(name='correct_count', initializer='zeros')
        self.total_count = self.add_weight(name='total_count', initializer='zeros')

    def update_state(self, y_true, y_pred):
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # remove the predictions and ground truth for the non-background classes
        y_pred = y_pred[:,:,:,0]
        y_true = y_true[:,:,:,0]
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.reduce_sum(tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32))
        total = tf.cast(tf.size(y_pred), dtype=tf.float32)
        self.correct_count.assign_add(correct)
        self.total_count.assign_add(total)

    def result(self):
        accuracy = (K.epsilon() + self.correct_count) / (K.epsilon() + self.total_count)#new
        return accuracy

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.correct_count.assign(0.0)
        self.total_count.assign(0.0)


# epochwise tissue accuracy
class TissueTypeAccuracy(tf.keras.metrics.Metric):
    def __init__(self, name='tissue_type_accuracy', **kwargs):
        super(TissueTypeAccuracy, self).__init__(name=name, **kwargs)
        self.correct_count = self.add_weight(name='correct_count', initializer='zeros')
        self.total_count = self.add_weight(name='total_count', initializer='zeros')

    def update_state(self, y_true, y_pred):
        y_true = tf.argmax(y_true, axis=-1)
        y_pred = tf.argmax(y_pred, axis=-1)
        
        tissue_indices = tf.cast(y_true, dtype=tf.bool)
        
        y_pred = tf.boolean_mask(y_pred, tissue_indices)
        y_true = tf.boolean_mask(y_true, tissue_indices)
        
        correct = tf.reduce_sum(tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32))
        total = tf.cast(tf.size(y_pred), dtype=tf.float32)
        self.correct_count.assign_add(correct)
        self.total_count.assign_add(total)

        
    def result(self):
        accuracy = (K.epsilon() + self.correct_count) / (self.total_count + K.epsilon())#new
        return accuracy

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.correct_count.assign(0.0)
        self.total_count.assign(0.0)


# precision, recall, and f1 for background below
class F1Macro0(tf.keras.metrics.Metric):
    def __init__(self, name='f1_macro_0', **kwargs):
        super(F1Macro0, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_positives = self.add_weight(name='fp', initializer='zeros')
        self.false_negatives = self.add_weight(name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * correct[:,:,:,0]))
        self.false_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * incorrect[:,:,:,0]))
        self.false_negatives.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,0]) * incorrect[:,:,:,0]))

    def result(self):
        precision = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_positives + K.epsilon())#new
        recall = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_negatives + K.epsilon())#new
        return 2 * (precision * recall) / (precision + recall)

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_positives.assign(0.0)
        self.false_negatives.assign(0.0)


class PrecisionMacro0(tf.keras.metrics.Metric):
    def __init__(self, name='precision_macro_0', **kwargs):
        super(PrecisionMacro0, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_positives = self.add_weight(name='fp', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * correct[:,:,:,0]))
        self.false_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * incorrect[:,:,:,0]))

    def result(self):
        precision = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_positives + K.epsilon())#new
        return precision

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_positives.assign(0.0)

class RecallMacro0(tf.keras.metrics.Metric):
    def __init__(self, name='recall_macro_0', **kwargs):
        super(RecallMacro0, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_negatives = self.add_weight(name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,0] * correct[:,:,:,0]))
        self.false_negatives.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,0]) * incorrect[:,:,:,0]))
        
    def result(self):
        recall = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_negatives + K.epsilon())#new
        return recall

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_negatives.assign(0,0)

# precision, recall, and f1 for neural below
class F1Macro1(tf.keras.metrics.Metric):
    def __init__(self, name='f1_macro_1', **kwargs):
        super(F1Macro1, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_positives = self.add_weight(name='fp', initializer='zeros')
        self.false_negatives = self.add_weight(name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * correct[:,:,:,1]))
        self.false_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * incorrect[:,:,:,1]))
        self.false_negatives.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,1]) * incorrect[:,:,:,1]))

    def result(self):
        precision = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_positives + K.epsilon())#new
        recall = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_negatives + K.epsilon())#new
        return 2 * (precision * recall) / (precision + recall)

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_positives.assign(0.0)
        self.false_negatives.assign(0.0)


class PrecisionMacro1(tf.keras.metrics.Metric):
    def __init__(self, name='precision_macro_1', **kwargs):
        super(PrecisionMacro1, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_positives = self.add_weight(name='fp', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * correct[:,:,:,1]))
        self.false_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * incorrect[:,:,:,1]))

    def result(self):
        precision = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_positives + K.epsilon())#new
        return precision

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_positives.assign(0.0)

class RecallMacro1(tf.keras.metrics.Metric):
    def __init__(self, name='recall_macro_1', **kwargs):
        super(RecallMacro1, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_negatives = self.add_weight(name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,1] * correct[:,:,:,1]))
        self.false_negatives.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,1]) * incorrect[:,:,:,1]))

    def result(self):
        recall = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_negatives + K.epsilon())#new
        return recall

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_negatives.assign(0.0)

# precision, recall, and f1 for mesenchyme below
class F1Macro2(tf.keras.metrics.Metric):
    def __init__(self, name='f1_macro_2', **kwargs):
        super(F1Macro2, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_positives = self.add_weight(name='fp', initializer='zeros')
        self.false_negatives = self.add_weight(name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * correct[:,:,:,2]))
        self.false_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * incorrect[:,:,:,2]))
        self.false_negatives.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,2]) * incorrect[:,:,:,2]))

    def result(self):
        precision = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_positives + K.epsilon())#new
        recall = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_negatives + K.epsilon())#new
        return 2 * (precision * recall) / (precision + recall)

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_positives.assign(0.0)
        self.false_negatives.assign(0.0)


class PrecisionMacro2(tf.keras.metrics.Metric):
    def __init__(self, name='precision_macro_2', **kwargs):
        super(PrecisionMacro2, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_positives = self.add_weight(name='fp', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * correct[:,:,:,2]))
        self.false_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * incorrect[:,:,:,2]))

    def result(self):
        precision = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_positives + K.epsilon())#new
        return precision

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_positives.assign(0.0)

class RecallMacro2(tf.keras.metrics.Metric):
    def __init__(self, name='recall_macro_2', **kwargs):
        super(RecallMacro2, self).__init__(name=name, **kwargs)
        self.true_positives = self.add_weight(name='tp', initializer='zeros')
        self.false_negatives = self.add_weight(name='fn', initializer='zeros')

    def update_state(self, y_true, y_pred):
        # convert y_pred from a matrix of softmax probabilities to one-hot encoded predictions
        y_pred_flat = tf.math.argmax(y_pred, axis = -1)
        y_pred = tf.one_hot(y_pred_flat, depth = y_pred.shape[-1])
        # construct matrix represention of where predictions matched the ground truth
        correct = tf.cast(tf.math.equal(y_pred, y_true), dtype=tf.float32)
        # construct matrix represention of where predictions didn't match the ground truth
        incorrect = 1.0 - correct
        # compute the true and false positives and negatives
        # for background
        self.true_positives.assign_add(tf.reduce_sum(y_pred[:,:,:,2] * correct[:,:,:,2]))
        self.false_negatives.assign_add(tf.reduce_sum((1.0 - y_pred[:,:,:,2]) * incorrect[:,:,:,2]))

    def result(self):
        recall = (K.epsilon() + self.true_positives) / (self.true_positives + self.false_negatives + K.epsilon())#new
        return recall

    def reset_states(self):
        # The state of the metric will be reset at the start of each epoch.
        self.true_positives.assign(0.0)
        self.false_negatives.assign(0.0)


########### below here are metric implementations from the MATLAB code

# evaluates important metrics for finetuning manually, without using keras
#note: it's assumed that the predictions and masks are flattened to a vector
def evaluate_metrics(preds, masks):
    
    classes = np.unique(masks)
    num_classes = len(classes)

    # get the confusion matrix
    tps = np.zeros([num_classes])
    tns = np.zeros([num_classes])
    fps = np.zeros([num_classes])
    fns = np.zeros([num_classes])
    total_class = np.zeros([num_classes])
    
    for i in range(num_classes):
        i_class = classes[i]
        
        tp = (preds == i_class) & (masks == i_class)
        tn = (preds != i_class) & (masks != i_class)
        fp = (preds == i_class) & (masks != i_class)
        fn = (preds != i_class) & (masks == i_class)
        
        tps[i] = np.sum(tp) + K.epsilon()
        tns[i] = np.sum(tn) + K.epsilon()
        fps[i] = np.sum(fp) + K.epsilon()
        fns[i] = np.sum(fn) + K.epsilon()
        total_class[i] = np.sum(masks == i_class)

    # compute the macroaveraged f1 score
    recalls = np.zeros([num_classes])
    precisions = np.zeros([num_classes])

    precisions = tps / (tps + fps)
    recalls = tps / (tps + fns)

    precision_M = np.sum(precisions) / num_classes
    recall_M = np.sum(recalls) / num_classes

    f1_score = 2 * ( (precision_M * recall_M) / (precision_M + recall_M) )
    
    # compute the binary accuracy
    binary_accuracy = (tps[0] + tns[0]) / (tps[0] + tns[0] + fps[0] + fns[0])

    # compute the tissue accuracy
    tissue_preds = preds[masks.astype(np.bool)]
    tissue_masks = masks[masks.astype(np.bool)]

    correct = np.sum(tissue_preds == tissue_masks)
    total = tissue_masks.size
    
    tissue_accuracy = correct / total

    return f1_score, binary_accuracy, tissue_accuracy
