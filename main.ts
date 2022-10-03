export function getPartialText(text: string): string {
 return text.split('').filter((_character, index) => index % 2 === 0).join('')
}
