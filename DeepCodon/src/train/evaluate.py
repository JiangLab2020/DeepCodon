import torch


def evaluate(model, data_loader, criterion, epoch):
    total_correct, total_loss, total_count = 0, 0, 0
    model.eval()
    with torch.no_grad():
        for enc_inputs, dec_inputs, dec_outputs in data_loader:
            enc_inputs, dec_inputs, dec_outputs = (
                enc_inputs.cuda(),
                dec_inputs.cuda(),
                dec_outputs.cuda(),
            )
            outputs = model(enc_inputs, dec_inputs)
            loss = criterion(outputs, dec_outputs.view(-1))
            total_loss += loss.item()
            non_zero_mask = dec_outputs.view(-1) != 0
            _, predicted = torch.max(outputs, 1)
            correct = (
                (predicted[non_zero_mask] == dec_outputs.view(-1)[non_zero_mask])
                .sum()
                .item()
            )
            total_correct += correct
            total_count += non_zero_mask.sum().item()
    avg_acc = total_correct / total_count
    avg_loss = total_loss / len(data_loader)
    return avg_loss, avg_acc
